window.addEventListener("DOMContentLoaded", function() {
  function normalizePath(path) {
    var normalized = [];
    path.split("/").forEach(function(bit, i) {
      if (bit === "." || (bit === "" && i !== 0)) {
        return;
      } else if (bit === "..") {
        if (normalized.length === 1 && normalized[0] === "") {
          // We must be trying to .. past the root!
          throw new Error("invalid path");
        } else if (normalized.length === 0 ||
                   normalized[normalized.length - 1] === "..") {
          normalized.push("..");
        } else {
          normalized.pop();
        }
      } else {
        normalized.push(bit);
      }
    });
    return normalized.join("/");
  }

  function insertAfter(referenceNode, newNode) {
    referenceNode.parentNode.insertBefore(newNode, referenceNode.nextSibling);
  }

  // `base_url` comes from the base.html template for this theme.
  var REL_BASE_URL = base_url;
  var ABS_BASE_URL = normalizePath(window.location.pathname + "/" +
                                   REL_BASE_URL);
  var CURRENT_VERSION = ABS_BASE_URL.split("/").pop();

  function makeSelect(options, selected) {
    var select = document.createElement("select");
    select.classList.add("form-control");

    options.forEach(function(i) {
      var option = new Option(i.text, i.value, undefined,
                              i.value === selected);
      select.add(option);
    });

    return select;
  }

  var xhr = new XMLHttpRequest();
  xhr.open("GET", REL_BASE_URL + "/../versions.json");
  xhr.onload = function() {

    try {
      var versions = JSON.parse(this.responseText);

      var realVersionDict = versions.find(function(i) {
        return i.version === CURRENT_VERSION ||
               i.aliases.includes(CURRENT_VERSION);
      });

      var realVersion = realVersionDict.version;

    } catch (error) {

      console.error(error);
      // Note - this branch is for local development server
      var realVersionDict = {"version": "local", "title": "local", "aliases": ["dev"]}
      var versions = [realVersionDict]
      var realVersion = "local"

    }

    if ( !realVersionDict.aliases.includes('latest') && !realVersionDict.aliases.includes('dev') ) {
      var dep_message = document.createElement("div");
      dep_message.id = "dep-message";
      dep_message.innerHTML = "This is the documentation for an old version of PALM. The current version can be found <a href=/latest" +
                              window.location.pathname.replace(/^\/?[^\/]+/, "") + ">here</a>.";
      document.body.insertBefore(dep_message, document.body.firstChild);
    }

    var select = makeSelect(versions.map(function(i) {
      let format_text = i.title;
      if (i.aliases.length !== 0) {
        format_text += ' (' + i.aliases.join(', ') + ')'
      }
      return {text: format_text, value: i.version};
    }), realVersion);
    select.addEventListener("change", function(event) {
      window.location.href = REL_BASE_URL + "/../" + this.value +
          window.location.pathname.replace(/^\/?[^\/]+/, "");
    });

    var container = document.createElement("a");
    container.id = "version-selector-custom";
    container.appendChild(select);

    var title = document.querySelector("a.navbar-brand");
    var height = window.getComputedStyle(title).getPropertyValue("height");
    container.style.height = height;
    container.style.marginRight = '16px';

    // Get the reference element
    var brand = document.getElementsByClassName("navbar-brand")[0]
    // brand.style.float = 'left';

    org_html = brand.outerHTML;
    new_html = "<div id='navbar-brand-div'>" + org_html + "</div>";
    brand.outerHTML = new_html;
    var brand_div = document.getElementById("navbar-brand-div")
    brand_div.appendChild(container);
    // Insert the new element into before sp2
    // parentDiv.insertBefore(sp1, sp2)
    // brand.parentNode.insertBefore(container, brand.nextSibling)
  };
  xhr.send();
});
