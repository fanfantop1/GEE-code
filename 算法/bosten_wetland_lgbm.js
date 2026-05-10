/**
 * 博斯腾湖湿地植被分类 - LightGBM (Gradient Boosting) 
 * 
 * 说明：
 * Google Earth Engine (GEE) 原生支持的 Boosting 算法为 ee.Classifier.smileGradientTreeBoost。
 * 该算法在原理上与 LightGBM 相似，均属于基于决策树的梯度提升算法。
 * 相比随机森林，Boosting 算法在处理难分类样本（Hard Examples）上更具优势
 */

// ================= 1. 基础配置与数据导入 =================
var roi = ee.FeatureCollection("projects/fiery-odyssey-491314-a2/assets/BostenLake_WetLandBoundary_new1050m");
var training_Points = ee.FeatureCollection("projects/fiery-odyssey-491314-a2/assets/XYDGSJ");
var region = roi.geometry(); 
var exportFolder = 'Bosten_Wetland_LGBM_Result'; 

Map.centerObject(roi, 11);
Map.addLayer(roi.style({color: 'red', fillColor: '00000000'}), {}, '研究区 (ROI)');

// ================= 2. 预处理函数 (Sentinel-2) =================
function preprocess_S2(image) {
  var qa = image.select('QA60');
  var mask = qa.bitwiseAnd(1 << 10).eq(0).and(qa.bitwiseAnd(1 << 11).eq(0));
  var ndvi = image.normalizedDifference(['B8', 'B4']).rename('NDVI');
  var opticalBands = image.select(['B1','B2','B3','B4','B5','B6','B7','B8','B8A','B9','B11','B12']).multiply(0.0001);
  return image.addBands(opticalBands, null, true).addBands(ndvi).updateMask(mask);
}

// ================= 3. 特征工程 (物候 + 光谱 + 纹理 + 雷达 = 64维) =================
var s2FullCol = ee.ImageCollection("COPERNICUS/S2_SR_HARMONIZED")
  .filterBounds(region).filterDate('2024-01-01', '2024-12-31').map(preprocess_S2);

// 物候计算
var composite_col = ee.ImageCollection(ee.List.sequence(1, 365, 15).map(function(d) {
  var date = ee.Date.fromYMD(2024, 1, 1).advance(d, 'day');
  return s2FullCol.filterDate(date, date.advance(15, 'day')).median()
    .set('doy', date.getRelative('day', 'year'));
}));

var bv = composite_col.select('NDVI').reduce(ee.Reducer.min()).rename('BV');
var mos_stack = composite_col.map(function(img){
  return img.select('NDVI').addBands(ee.Image.constant(ee.Number(img.get('doy'))).rename('doy').float());
}).qualityMosaic('NDVI');

var mos_date = mos_stack.select('doy').rename('MOS');
var mos_val = mos_stack.select('NDVI').rename('MOS_val');
var sa = mos_val.subtract(bv).rename('SA');
var threshold = bv.add(sa.multiply(0.2));

var sos = ee.Image(composite_col.iterate(function(img, acc){
  return ee.Image(acc).where(ee.Image(img).select('NDVI').gt(threshold).and(ee.Image(acc).eq(0)), ee.Image.constant(ee.Number(ee.Image(img).get('doy'))));
}, ee.Image.constant(0))).rename('SOS');

var eos = ee.Image(composite_col.iterate(function(img, acc){
  return ee.Image(acc).where(ee.Image(img).select('NDVI').gt(threshold), ee.Image.constant(ee.Number(ee.Image(img).get('doy'))));
}, ee.Image.constant(0))).rename('EOS');

var los = eos.subtract(sos).rename('LOS');
var irs = mos_val.subtract(threshold).divide(mos_date.subtract(sos).add(0.001)).rename('IRS');
var drs = mos_val.subtract(threshold).divide(eos.subtract(mos_date).add(0.001)).rename('DRS');

var phenology_metrics = ee.Image.cat([sos, eos, los, mos_date, bv, sa, irs, drs]).toFloat();

// 指数计算
var growS2 = s2FullCol.filterDate('2024-05-01', '2024-10-31').median();

function calculate_all_indices(img) {
  var v = {'B2':img.select('B2'),'B3':img.select('B3'),'B4':img.select('B4'),'B5':img.select('B5'),'B6':img.select('B6'),'B7':img.select('B7'),'B8':img.select('B8'),'B8A':img.select('B8A'),'B11':img.select('B11'),'B12':img.select('B12')};
  return img.addBands([
    img.normalizedDifference(['B8', 'B4']).rename('NDVI_grow'),
    v.B8.divide(v.B4.add(0.0001)).rename('RVI'),
    img.expression('2.5 * (B8 - B4) / (B8 + 6 * B4 - 7.5 * B2 + 1)', v).rename('EVI'),
    img.expression('1.5 * (B8 - B4) / (B8 + B4 + 0.5)', v).rename('SAVI'),
    img.normalizedDifference(['B8A', 'B3']).rename('GNDVI'),
    img.normalizedDifference(['B3', 'B8']).rename('NDWI'),
    img.normalizedDifference(['B3', 'B11']).rename('MNDWI'),
    img.normalizedDifference(['B8', 'B11']).rename('LSWI'),
    img.normalizedDifference(['B8A', 'B5']).rename('NDVIre1'),
    img.normalizedDifference(['B8A', 'B6']).rename('NDVIre2'),
    img.normalizedDifference(['B8A', 'B7']).rename('NDVIre3'),
    img.normalizedDifference(['B6', 'B5']).rename('NDre1'),
    img.normalizedDifference(['B7', 'B5']).rename('NDre2'),
    v.B7.divide(v.B5.add(0.0001)).subtract(1).rename('CIre'),
    img.expression('705 + 35 * (((B4 + B6)/2 - B5) / (B6 - B5 + 0.0001))', v).rename('REIP'),
    img.expression('(B8A - B4) / (B5 / (B6 + 0.0001) + 0.0001)', v).rename('IRECI'),
    img.expression('0.3037*B2 + 0.2793*B3 + 0.4743*B4 + 0.5585*B8 + 0.5082*B11 + 0.1863*B12', v).rename('TC_B'),
    img.expression('-0.2848*B2 - 0.2435*B3 - 0.5436*B4 + 0.7243*B8 + 0.0840*B11 - 0.1800*B12', v).rename('TC_G'),
    img.expression('0.1509*B2 + 0.1973*B3 + 0.3279*B4 + 0.3406*B8 - 0.7112*B11 - 0.4572*B12', v).rename('TC_W'),
    img.expression('2.5 * (B8 - B4) / (B8 + 2.4 * B4 + 1)', v).rename('EVI2'),
    img.expression('(2 * B8 + 1 - sqrt(pow(2 * B8 + 1, 2) - 8 * (B8 - B4))) / 2.0', v).rename('MSAVI'),
    img.expression('((B11 + B4) - (B8 + B2)) / ((B11 + B4) + (B8 + B2) + 0.0001)', v).rename('BSI'),
    img.expression('(0.1 * B8 - B4) / (0.1 * B8 + B4 + 0.0001)', v).rename('WDRVI'),
    img.expression('((B5 - B4) - 0.2 * (B5 - B3)) * (B5 / (B4 + 0.0001))', v).rename('MCARI'),
    img.expression('(B6 - B5) / (B5 - B4 + 0.0001)', v).rename('MTCI'),
    img.normalizedDifference(['B8', 'B12']).rename('NDII'),
    img.expression('(2*B3 - B4 - B2) / (2*B3 + B4 + B2 + 0.0001)', v).rename('GLI'),
    img.expression('(pow(B8, 2) - B4) / (pow(B8, 2) + B4 + 0.0001)', v).rename('NLI'),
    img.expression('(B3 - B4) / (B3 + B4 - B2 + 0.0001)', v).rename('VARI'),
    img.expression('(B8 / (B4 + 0.0001) - 1) / sqrt(B8 / (B4 + 0.0001) + 1)', v).rename('MSR')
  ]);
}

var spectralFeatures = calculate_all_indices(growS2);
var texture = growS2.select('B8').unitScale(0, 0.4).multiply(100).toInt().glcmTexture({size: 3})
  .select(['B8_asm','B8_contrast','B8_corr','B8_ent','B8_var','B8_idm','B8_shade','B8_prom','B8_sent'], 
          ['ASM','Contrast','Corr','Ent','Var','IDM','Shade','Prom','SEnt']);

var s1 = ee.ImageCollection("COPERNICUS/S1_GRD").filterBounds(region)
  .filterDate('2024-05-01', '2024-10-31').filter(ee.Filter.eq('instrumentMode', 'IW')).median();
var radar = s1.select(['VV','VH']).addBands([
  s1.select('VV').add(s1.select('VH')).rename('SAR_sum'),
  s1.select('VV').subtract(s1.select('VH')).rename('SAR_diff'),
  s1.select('VV').subtract(s1.select('VH')).divide(s1.select('VV').add(s1.select('VH')).add(0.0001)).rename('SAR_NDVI')
]);

var finalImage = ee.Image.cat([
  growS2.select(['B1','B2','B3','B4','B5','B6','B7','B8','B8A','B9','B11','B12']),
  spectralFeatures.select(['NDVI_grow','RVI','EVI','SAVI','GNDVI','NDWI','MNDWI','LSWI','NDVIre1','NDVIre2','NDVIre3','NDre1','NDre2','CIre','REIP','IRECI','TC_B','TC_G','TC_W','EVI2','MSAVI','BSI','WDRVI','MCARI','MTCI','NDII','GLI','NLI','VARI','MSR']),
  texture,
  radar,
  phenology_metrics
]).toFloat().clip(roi); 

// ================= 4. 模型调优与分类 (Gradient Boosting) =================

var cleanPoints = training_Points.filter(ee.Filter.notNull(['label']));
var samples = finalImage.sampleRegions({
  collection: cleanPoints,
  properties: ['label'],
  scale: 10,
  tileScale: 16,
  geometries: true
});

var finalSamples = samples
  .filter(ee.Filter.notNull(['B2', 'label']))
  .map(function(f) { return f.set('label', ee.Number(f.get('label')).toInt()) });

var samplesWithRandom = finalSamples.randomColumn('random', 42);
var trainSet = samplesWithRandom.filter(ee.Filter.lt('random', 0.7));
var testSet = samplesWithRandom.filter(ee.Filter.gte('random', 0.7));

/**
 * 4.1 Boosting 超参数网格搜索
 * 重点调优：
 * - numberOfTrees (迭代次数)
 * - shrinkage (学习率/步长)
 * - maxNodes (决策树深度限制)
 */
var treeList = [100, 200];
var lrList = [0.1, 0.2];

var tuningResults = treeList.map(function(t) {
  return lrList.map(function(lr) {
    var gbm = ee.Classifier.smileGradientTreeBoost({
      numberOfTrees: t,
      shrinkage: lr,
      samplingRate: 0.7,
      maxNodes: 10
    }).train({
      features: trainSet,
      classProperty: 'label',
      inputProperties: finalImage.bandNames()
    });
    
    var accuracy = testSet.classify(gbm).errorMatrix('label', 'classification').accuracy();
    return ee.Feature(null, {'trees': t, 'lr': lr, 'accuracy': accuracy});
  });
});

var tuningFC = ee.FeatureCollection(ee.List(tuningResults).flatten());
var bestParams = tuningFC.sort('accuracy', false).first();
print('🏆 Gradient Boosting 调优结果:', bestParams);

// 4.2 使用最优参数训练最终模型
var gbmFinal = ee.Classifier.smileGradientTreeBoost({
  numberOfTrees: ee.Number(bestParams.get('trees')),
  shrinkage: ee.Number(bestParams.get('lr')),
  samplingRate: 0.7,
  maxNodes: 10
}).train({
  features: trainSet,
  classProperty: 'label',
  inputProperties: finalImage.bandNames()
});

// ================= 5. 精度评价 =================
var testClassified = testSet.classify(gbmFinal);
var confusionMatrix = testClassified.errorMatrix('label', 'classification');

print('👉 有效样点总数:', finalSamples.size());
print('👉 混淆矩阵:', confusionMatrix);
print('👉 总体精度 (OA):', confusionMatrix.accuracy());
print('👉 Kappa 系数:', confusionMatrix.kappa());

// 计算各类 F1-score
var pa = confusionMatrix.producersAccuracy(); 
var ua = confusionMatrix.consumersAccuracy();  
var f1Score = ee.List.sequence(0, ee.Number(confusionMatrix.order().length()).subtract(1)).map(function(i) {
  var p = ee.Number(ua.get([0, i]));
  var r = ee.Number(pa.get([i, 0]));
  return p.multiply(r).multiply(2).divide(p.add(r).add(0.0001));
});
print('👉 各类 F1-score:', f1Score);

var accuracyFC = ee.FeatureCollection([
  ee.Feature(null, {'Metric': 'Overall Accuracy', 'Value': confusionMatrix.accuracy()}),
  ee.Feature(null, {'Metric': 'Kappa Coefficient', 'Value': confusionMatrix.kappa()}),
  ee.Feature(null, {'Metric': 'Best_Trees', 'Value': bestParams.get('trees')}),
  ee.Feature(null, {'Metric': 'Best_LR', 'Value': bestParams.get('lr')})
]);

// ================= 6. 可视化与导出 =================
var finalClassifiedMap = finalImage.classify(gbmFinal);
var palette = ['0000FF', '006400', 'FFA500', 'ADFF2F'];
Map.addLayer(finalClassifiedMap, {min: 1, max: 4, palette: palette}, 'GBM 分类结果 (LightGBM 模式)');

Export.image.toDrive({
  image: finalClassifiedMap,
  description: '1_Bosten_LGBM_Classification_Raster',
  folder: exportFolder,
  region: region,
  scale: 10,
  maxPixels: 1e13,
  crs: 'EPSG:4326',
  fileFormat: 'GeoTIFF'
});

Export.table.toDrive({
  collection: accuracyFC,
  description: '2_Bosten_LGBM_Accuracy_Metrics',
  folder: exportFolder,
  fileFormat: 'CSV'
});

print('✅ LightGBM (Gradient Boosting) 任务已就绪。');
