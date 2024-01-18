document.addEventListener('DOMContentLoaded', (event) => {
  document.querySelectorAll('footer p').forEach((el) => {
    el.parentNode.removeChild(el);
  });
});
