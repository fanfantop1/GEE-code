/**
 * 2024年夏季最佳哨兵1/哨兵2影像下载脚本
 * =====================================
 * 输出：
 * 1. Sentinel-2 真彩色影像（RGB = B4/B3/B2），选择最佳单日并对当日全部影像拼接
 * 2. Sentinel-1 夏季 VV/VH 双极化影像，选择最佳单日并对当日全部影像拼接
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
var targetDate = ee.Date.fromYMD(year, 7, 15);
var driveFolder = 'Best_S1_S2_2024_Summer';
var crs = 'EPSG:4326';
var scaleS2 = 10;
var scaleS1 = 10;

Map.centerObject(roi, 10);
Map.addLayer(roi, {color: 'red'}, 'ROI');

// ===================== 通用函数 =====================
function addDateString(image) {
  return image.set('date_str', image.date().format('YYYY-MM-dd'));
}

function getMaskCoverage(image, bandName, region, scale) {
  var validMask = image.select(bandName).mask();
  var validCount = ee.Number(validMask.reduceRegion({
    reducer: ee.Reducer.sum(),
    geometry: region,
    scale: scale,
    maxPixels: 1e9
  }).get(bandName));

  var totalCount = ee.Number(ee.Image.constant(1).clip(region).reduceRegion({
    reducer: ee.Reducer.count(),
    geometry: region,
    scale: scale,
    maxPixels: 1e9
  }).get('constant'));

  return ee.Number(ee.Algorithms.If(
    totalCount.gt(0),
    validCount.divide(totalCount),
    0
  ));
}

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

var s2Raw = ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED')
  .filterBounds(roi)
  .filterDate(startDate, endDate)
  .filter(ee.Filter.lte('CLOUDY_PIXEL_PERCENTAGE', 20))
  .map(addDateString);

var s2Dates = ee.List(s2Raw.aggregate_array('date_str')).distinct().sort();
var s2Daily = ee.ImageCollection.fromImages(s2Dates.map(function(dateStr) {
  dateStr = ee.String(dateStr);
  var dayStart = ee.Date.parse('YYYY-MM-dd', dateStr);
  var dayEnd = dayStart.advance(1, 'day');
  var dayCol = s2Raw.filterDate(dayStart, dayEnd);
  var mosaic = dayCol.map(maskS2).mosaic().clip(roi);
  var validRatio = getMaskCoverage(mosaic, 'B4', roi, 20);
  var meanCloud = ee.Number(dayCol.aggregate_mean('CLOUDY_PIXEL_PERCENTAGE'));
  var sceneCount = dayCol.size();
  var daysFromTarget = ee.Number(dayStart.difference(targetDate, 'day')).abs();
  var score = validRatio.multiply(1000)
    .subtract(meanCloud.multiply(2))
    .subtract(daysFromTarget.multiply(0.1));

  return mosaic.set({
    date_str: dateStr,
    valid_ratio: validRatio,
    mean_cloud_percent: meanCloud,
    scene_count: sceneCount,
    days_from_target: daysFromTarget,
    quality_score: score
  });
})).sort('quality_score', false);

var bestS2 = ee.Image(s2Daily.first());

print('Sentinel-2 候选原始影像数量:', s2Raw.size());
print('Sentinel-2 候选日期数量:', s2Dates.size());
print('最佳 Sentinel-2 单日拼接信息:', bestS2);
print('最佳 Sentinel-2 日期:', bestS2.get('date_str'));
print('最佳 Sentinel-2 平均云量(%):', bestS2.get('mean_cloud_percent'));
print('最佳 Sentinel-2 有效像元比例:', bestS2.get('valid_ratio'));
print('最佳 Sentinel-2 当日影像数:', bestS2.get('scene_count'));

Map.addLayer(
  bestS2,
  {bands: ['B4', 'B3', 'B2'], min: 0, max: 3000},
  'Best Sentinel-2 RGB'
);

// ===================== Sentinel-1 =====================
// 哨兵1无云影响，这里选择最佳单日并对该日所有影像拼接：
// 1. 优先研究区覆盖更完整
// 2. 覆盖接近时优先靠近夏季中期
var s1Raw = ee.ImageCollection('COPERNICUS/S1_GRD')
  .filterBounds(roi)
  .filterDate(startDate, endDate)
  .filter(ee.Filter.eq('instrumentMode', 'IW'))
  .filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VV'))
  .filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VH'))
  .select(['VV', 'VH'])
  .map(addDateString);

var s1Dates = ee.List(s1Raw.aggregate_array('date_str')).distinct().sort();
var s1Daily = ee.ImageCollection.fromImages(s1Dates.map(function(dateStr) {
  dateStr = ee.String(dateStr);
  var dayStart = ee.Date.parse('YYYY-MM-dd', dateStr);
  var dayEnd = dayStart.advance(1, 'day');
  var dayCol = s1Raw.filterDate(dayStart, dayEnd);
  var mosaic = dayCol.mosaic().clip(roi);
  var coverageRatio = getMaskCoverage(mosaic, 'VV', roi, 20);
  var sceneCount = dayCol.size();
  var daysFromTarget = ee.Number(dayStart.difference(targetDate, 'day')).abs();
  var score = coverageRatio.multiply(1000)
    .subtract(daysFromTarget.multiply(0.1));

  return mosaic.set({
    date_str: dateStr,
    coverage_ratio: coverageRatio,
    scene_count: sceneCount,
    days_from_target: daysFromTarget,
    quality_score: score,
    orbit_pass_list: ee.List(dayCol.aggregate_array('orbitProperties_pass')).distinct().join('|'),
    relative_orbit_list: ee.List(dayCol.aggregate_array('relativeOrbitNumber_start')).distinct().join('|')
  });
})).sort('quality_score', false);

var bestS1 = ee.Image(s1Daily.first());

print('Sentinel-1 候选原始影像数量:', s1Raw.size());
print('Sentinel-1 候选日期数量:', s1Dates.size());
print('最佳 Sentinel-1 单日拼接信息:', bestS1);
print('最佳 Sentinel-1 日期:', bestS1.get('date_str'));
print('最佳 Sentinel-1 覆盖比例:', bestS1.get('coverage_ratio'));
print('最佳 Sentinel-1 当日影像数:', bestS1.get('scene_count'));
print('最佳 Sentinel-1 距离目标日期天数:', bestS1.get('days_from_target'));

// ===================== 影像信息汇总输出 =====================
var imageInfo = ee.FeatureCollection([
  ee.Feature(null, {
    sensor: 'Sentinel-2',
    mosaic_type: 'best_single_day_mosaic',
    date: bestS2.get('date_str'),
    cloud_percent: bestS2.get('mean_cloud_percent'),
    valid_ratio: bestS2.get('valid_ratio'),
    scene_count: bestS2.get('scene_count'),
    days_from_target: bestS2.get('days_from_target')
  }),
  ee.Feature(null, {
    sensor: 'Sentinel-1',
    mosaic_type: 'best_single_day_mosaic',
    date: bestS1.get('date_str'),
    coverage_ratio: bestS1.get('coverage_ratio'),
    scene_count: bestS1.get('scene_count'),
    days_from_target: bestS1.get('days_from_target'),
    orbit_pass: bestS1.get('orbit_pass_list'),
    relative_orbit: bestS1.get('relative_orbit_list')
  })
]);

print('最佳影像信息汇总:', imageInfo);

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

Export.table.toDrive({
  collection: imageInfo,
  description: 'Best_S1_S2_2024_Summer_Info',
  folder: driveFolder,
  fileNamePrefix: 'Best_S1_S2_2024_Summer_Info',
  fileFormat: 'CSV'
});

print('已创建 3 个导出任务：Sentinel-2 真彩色 + Sentinel-1 VV/VH + 影像信息汇总 CSV。');
print('请到 GEE 右侧 Tasks 面板点击 Run。');
