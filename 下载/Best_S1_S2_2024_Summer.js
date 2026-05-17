/**
 * 2024年夏季最佳哨兵1/哨兵2影像下载脚本
 * =====================================
 * 输出：
 * 1. Sentinel-2 真彩色影像（RGB = B4/B3/B2），优先选择研究区内有效像元比例高、云量低的单景
 * 2. Sentinel-1 夏季 VV/VH 双极化影像，选择单张最佳影像而不是时间合成
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
// 哨兵1无云影响，这里选择夏季中的最佳单景：
// 1. 优先研究区覆盖更完整的影像
// 2. 在覆盖接近时，优先更靠近夏季中期（7月15日）的影像
function addS1Quality(image) {
  var footprint = ee.Geometry(image.geometry());
  var roiArea = roi.area(1);
  var intersectArea = footprint.intersection(roi, 1).area(1);
  var coverageRatio = ee.Algorithms.If(
    roiArea.gt(0),
    intersectArea.divide(roiArea),
    0
  );

  var daysFromTarget = ee.Number(image.date().difference(targetDate, 'day')).abs();

  return image
    .set('coverage_ratio', coverageRatio)
    .set('days_from_target', daysFromTarget);
}

var s1Col = ee.ImageCollection('COPERNICUS/S1_GRD')
  .filterBounds(roi)
  .filterDate(startDate, endDate)
  .filter(ee.Filter.eq('instrumentMode', 'IW'))
  .filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VV'))
  .filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VH'))
  .select(['VV', 'VH'])
  .map(addS1Quality)
  .sort('days_from_target')
  .sort('coverage_ratio', false);

var bestS1 = ee.Image(s1Col.first()).clip(roi);

print('Sentinel-1 候选影像数量:', s1Col.size());
print('最佳 Sentinel-1 影像信息:', bestS1);
print('最佳 Sentinel-1 日期:', ee.Date(bestS1.get('system:time_start')).format('YYYY-MM-dd'));
print('最佳 Sentinel-1 覆盖比例:', bestS1.get('coverage_ratio'));
print('最佳 Sentinel-1 距离目标日期天数:', bestS1.get('days_from_target'));

// ===================== 影像信息汇总输出 =====================
var imageInfo = ee.FeatureCollection([
  ee.Feature(null, {
    sensor: 'Sentinel-2',
    image_id: bestS2.get('PRODUCT_ID'),
    date: ee.Date(bestS2.get('system:time_start')).format('YYYY-MM-dd'),
    cloud_percent: bestS2.get('CLOUDY_PIXEL_PERCENTAGE'),
    valid_ratio: bestS2.get('valid_ratio'),
    orbit: bestS2.get('SENSING_ORBIT_NUMBER'),
    tile: bestS2.get('MGRS_TILE')
  }),
  ee.Feature(null, {
    sensor: 'Sentinel-1',
    image_id: bestS1.get('system:index'),
    date: ee.Date(bestS1.get('system:time_start')).format('YYYY-MM-dd'),
    coverage_ratio: bestS1.get('coverage_ratio'),
    days_from_target: bestS1.get('days_from_target'),
    orbit_pass: bestS1.get('orbitProperties_pass'),
    relative_orbit: bestS1.get('relativeOrbitNumber_start')
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
