/**
 * 2024年夏季最佳哨兵1/哨兵2影像下载脚本
 * =====================================
 * 输出：
 * 1. Sentinel-2 真彩色影像（RGB = B4/B3/B2），优先选择研究区内有效像元比例高、云量低的影像
 * 2. Sentinel-1 夏季 VV/VH 双极化影像（使用夏季中值合成得到一张高质量影像）
 *
 * 使用前请修改：
 * 1. roi - 研究区，可替换为你自己的 Geometry 或 Asset
 * 2. driveFolder - Google Drive 导出文件夹名称
 */

// ===================== 研究区配置 =====================
// 示例区域：博斯腾湖附近，可按需替换
var roi = ee.Geometry.Rectangle([86.5, 41.7, 87.5, 42.2]);

// 如果你有自己的矢量资产，可改成下面这种写法：
// var roi = ee.FeatureCollection('users/yourname/your_roi').geometry();

// ===================== 参数配置 =====================
var year = 2024;
var startDate = ee.Date.fromYMD(year, 6, 1);
var endDate = ee.Date.fromYMD(year, 9, 1);
var driveFolder = 'Best_S1_S2_2024_Summer';
var crs = 'EPSG:4326';
var scaleS2 = 10;
var scaleS1 = 10;

Map.centerObject(roi, 10);
Map.addLayer(roi, {color: 'red'}, 'ROI');

// ===================== Sentinel-2 =====================
// 使用 SCL 去除云、云影、卷云、雪冰，保留高质量像元
function maskS2(image) {
  var scl = image.select('SCL');
  var mask = scl.eq(4)   // vegetation
    .or(scl.eq(5))       // bare soil
    .or(scl.eq(6))       // water
    .or(scl.eq(7))       // unclassified
    .or(scl.eq(11));     // snow/ice, 如不需要可删除

  return image
    .updateMask(mask)
    .select(['B2', 'B3', 'B4', 'B8', 'SCL'])
    .copyProperties(image, image.propertyNames());
}

// 计算研究区内有效像元比例，用于挑选“最佳”影像
function addValidCoverage(image) {
  var validMask = image.select('B4').mask();
  var validCount = ee.Number(validMask.reduceRegion({
    reducer: ee.Reducer.sum(),
    geometry: roi,
    scale: 20,
    maxPixels: 1e9
  }).get('B4'));

  var totalCount = ee.Number(ee.Image.constant(1).clip(roi).reduceRegion({
    reducer: ee.Reducer.count(),
    geometry: roi,
    scale: 20,
    maxPixels: 1e9
  }).get('constant'));

  var validRatio = ee.Algorithms.If(
    totalCount.gt(0),
    validCount.divide(totalCount),
    0
  );

  return image.set('valid_ratio', validRatio);
}

var s2Col = ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED')
  .filterBounds(roi)
  .filterDate(startDate, endDate)
  .filter(ee.Filter.lte('CLOUDY_PIXEL_PERCENTAGE', 20))
  .map(maskS2)
  .map(addValidCoverage)
  .sort('valid_ratio', false)
  .sort('CLOUDY_PIXEL_PERCENTAGE');

var bestS2 = ee.Image(s2Col.first()).clip(roi);

print('Sentinel-2 候选影像数量:', s2Col.size());
print('最佳 Sentinel-2 影像信息:', bestS2);
print('最佳 Sentinel-2 日期:', ee.Date(bestS2.get('system:time_start')).format('YYYY-MM-dd'));
print('最佳 Sentinel-2 云量(%):', bestS2.get('CLOUDY_PIXEL_PERCENTAGE'));
print('最佳 Sentinel-2 有效像元比例:', bestS2.get('valid_ratio'));

Map.addLayer(
  bestS2,
  {bands: ['B4', 'B3', 'B2'], min: 0, max: 3000},
  'Best Sentinel-2 RGB'
);

// ===================== Sentinel-1 =====================
// 哨兵1无云影响，使用夏季 VV/VH 中值合成得到一张稳定高质量影像
var s1Col = ee.ImageCollection('COPERNICUS/S1_GRD')
  .filterBounds(roi)
  .filterDate(startDate, endDate)
  .filter(ee.Filter.eq('instrumentMode', 'IW'))
  .filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VV'))
  .filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VH'))
  .select(['VV', 'VH']);

var bestS1 = s1Col.median().clip(roi);

print('Sentinel-1 候选影像数量:', s1Col.size());
print('Sentinel-1 夏季合成影像:', bestS1);

Map.addLayer(
  bestS1,
  {bands: ['VH', 'VV', 'VH'], min: -25, max: 5},
  'Sentinel-1 VH/VV'
);

// ===================== 导出任务 =====================
Export.image.toDrive({
  image: bestS2.select(['B4', 'B3', 'B2']),
  description: 'Best_Sentinel2_RGB_2024_Summer',
  folder: driveFolder,
  fileNamePrefix: 'Best_Sentinel2_RGB_2024_Summer',
  region: roi,
  scale: scaleS2,
  crs: crs,
  maxPixels: 1e13
});

Export.image.toDrive({
  image: bestS1.select(['VV', 'VH']),
  description: 'Best_Sentinel1_VV_VH_2024_Summer',
  folder: driveFolder,
  fileNamePrefix: 'Best_Sentinel1_VV_VH_2024_Summer',
  region: roi,
  scale: scaleS1,
  crs: crs,
  maxPixels: 1e13
});

print('已创建 2 个导出任务：Sentinel-2 真彩色 + Sentinel-1 VV/VH。');
print('请到 GEE 右侧 Tasks 面板点击 Run。');
