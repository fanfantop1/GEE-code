/**
 * 哨兵2 长时序 NDVI 自动下载脚本 (2015-2025)
 * =============================================
 * 输出结构：
 *   Google Drive
 *   ├── NDVI_Annual/    → 逐年 NDVI（年均值 + 年最大值），11 幅
 *   └── NDVI_Monthly/   → 逐月 NDVI（2015-01 到 2025-12），132 幅
 *
 * 注意：
 *   - 哨兵2A 2015年6月才发射，2015-01~2015-05 的文件会是空栅格，可忽略
 *   - NDVI 存储值 = (真实NDVI × 100) + 100，还原公式: (DN - 100) / 100
 *
 * 使用前请修改：
 *   1. roi — 你的研究区矢量（可手绘 Geometry 或导入 Asset）
 *   2. startYear / endYear — 起止年份
 */

// ===================== 研究区配置 =====================
// 方式一：手绘矩形（示例：博斯腾湖区域）
var roi = ee.Geometry.Rectangle([86.5, 41.7, 87.5, 42.2]);

// 方式二：从 Assets 导入矢量
// var roi = ee.FeatureCollection("users/yourname/bosten_lake");

// ===================== 参数配置 =====================
var startYear  = 2015;
var endYear    = 2025;
var scale      = 10;          // 哨兵2 原始分辨率 10m
var crs        = 'EPSG:4326'; // 导出投影

// Google Drive 根文件夹（两个子文件夹会自动创建在其下）
var driveFolder = 'Sentinel2_NDVI';

// ===================== 云掩膜函数 =====================
// 使用 SCL 波段去除云/云影/雪
function maskS2clouds(image) {
  var scl = image.select('SCL');
  // SCL 有效像元: 2-7, 11 (植被/土壤/水体/裸地等)
  var mask = scl.eq(2).or(scl.eq(3)).or(scl.eq(4))
                  .or(scl.eq(5)).or(scl.eq(6))
                  .or(scl.eq(7)).or(scl.eq(11));
  return image.updateMask(mask);
}

// NDVI 计算
function addNDVI(image) {
  var ndvi = image.expression(
    '(NIR - RED) / (NIR + RED)', {
      'NIR': image.select('B8'),
      'RED': image.select('B4')
  }).rename('NDVI');
  return image.addBands(ndvi);
}

// ===================== 构建哨兵2 NDVI 集合 =====================
// S2 L2A 地表反射率，Harmonized 整合了新旧哨兵2
var s2 = ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED')
  .filterBounds(roi)
  .filter(ee.Filter.calendarRange(startYear, endYear, 'year'))
  .map(maskS2clouds)
  .map(addNDVI);

print('哨兵2 影像总数:', s2.size());

// ===================== 逐年 NDVI 计算与导出 =====================
// 每年生成一幅影像，包含 avg 和 max 两个波段
for (var y = startYear; y <= endYear; y++) {
  var yearCol = s2.filter(ee.Filter.calendarRange(y, y, 'year'));

  var ndviAvg = yearCol.select('NDVI').mean().rename('NDVI_avg');
  var ndviMax = yearCol.select('NDVI').max().rename('NDVI_max');

  var yearImage = ndviAvg.addBands(ndviMax)
    .set('year', y)
    .clip(roi);

  // 转为 Byte 类型 (-1~1 缩放至 0~200) 以减小文件体积
  yearImage = yearImage.multiply(100).add(100).toByte();

  Export.image.toDrive({
    image: yearImage,
    description: 'NDVI_Annual_' + y,
    folder: driveFolder + '/NDVI_Annual',
    fileNamePrefix: 'NDVI_' + y,
    region: roi,
    scale: scale,
    crs: crs,
    maxPixels: 1e13
  });
}

print('已提交逐年 NDVI 导出: ' + (endYear - startYear + 1) + ' 幅');

// ===================== 逐月 NDVI 计算与导出 =====================
// 每年每月各出一幅，2015-01 → 2025-12 共 132 幅
var monthCount = 0;

for (var y = startYear; y <= endYear; y++) {
  for (var m = 1; m <= 12; m++) {
    var monthCol = s2.filter(ee.Filter.calendarRange(y, y, 'year'))
                     .filter(ee.Filter.calendarRange(m, m, 'month'));

    var padMonth = m < 10 ? '0' + m : '' + m;
    var label = y + '-' + padMonth;

    var ndviMonthly = monthCol.select('NDVI').mean()
      .rename('NDVI')
      .set('year', y)
      .set('month', m)
      .clip(roi);

    ndviMonthly = ndviMonthly.multiply(100).add(100).toByte();

    Export.image.toDrive({
      image: ndviMonthly,
      description: 'NDVI_' + label,
      folder: driveFolder + '/NDVI_Monthly',
      fileNamePrefix: 'NDVI_' + label,
      region: roi,
      scale: scale,
      crs: crs,
      maxPixels: 1e13
    });

    monthCount++;
  }
}

print('已提交逐月 NDVI 导出: ' + monthCount + ' 幅');
print('全部完成！总共 ' + (endYear - startYear + 1 + monthCount) + ' 个导出任务，请打开 Google Drive 查看。');
print('提示：2015年1-5月（哨兵2A发射前）的文件为空栅格，可忽略。');
