/**
 * GEE Tasks 面板自动点击 Run + OK 脚本
 * =====================================
 * 使用：F12 → Console，粘贴全部代码，回车
 * 停止：GEE_Stop() 或 刷新页面
 */

window.GEE_Stop = function () {
  clearInterval(window.__geeTimer);
  console.log('已停止');
};

(function () {
  var interval = 1500;

  // 尝试从 shadow DOM 里找元素
  function deepQuery(selector) {
    // 先普通查找
    var el = document.querySelector(selector);
    if (el) return el;

    // 穿透所有 shadowRoot
    function walk(root) {
      if (!root) return null;
      try { var found = root.querySelector(selector); if (found) return found; } catch (e) {}
      var children = root.querySelectorAll ? root.querySelectorAll('*') : [];
      for (var i = 0; i < children.length; i++) {
        var sr = children[i].shadowRoot;
        if (sr) { var f = walk(sr); if (f) return f; }
      }
      return null;
    }
    return walk(document.documentElement);
  }

  window.__geeTimer = setInterval(function () {
    var runBtn = deepQuery('ee-button.run-button');
    var okBtn  = deepQuery('ee-button.ok-button');

    if (okBtn) {
      okBtn.click();
      console.log('已点击 OK');
    } else if (runBtn) {
      runBtn.click();
      console.log('已点击 Run');
    } else {
      // 顺便排查一下页面状态
      var taskPane = document.querySelector('.task-pane, [class*="task-pane"]');
      if (!taskPane) {
        console.log('⚠ 找不到 Tasks 面板，请先点击右侧 Tasks 标签展开');
      } else {
        console.log('没有待处理的任务了（可能已全部提交或正在运行中）');
      }
      // 不自动停止，持续检测（新任务可能会陆续出现）
    }
  }, interval);

  console.log('已开始自动点击（间隔 ' + interval + 'ms），输入 GEE_Stop() 可停止');
})();
