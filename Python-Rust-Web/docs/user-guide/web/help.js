(() => {
  "use strict";

  const storageKey = "magnetictb.web.locale";
  const body = document.body;

  // Older bookmarks now open the dedicated page instead of a long home section.
  if (body.dataset.page === "home") {
    const movedSections = {
      "installation": "guide/Installation.html",
      "conda": "guide/Installation.html",
      "guide": "guide/FirstModel.html",
      "symham": "guide/FirstModel.html",
      "tutorials": "tutorial/HamiltonianOutputs.html",
      "graphene": "tutorial/HamiltonianOutputs.html",
      "simple-cubic": "tutorial/HamiltonianOutputs.html",
      "mos2": "tutorial/HamiltonianOutputs.html",
      "reference": "guide/FirstModel.html",
      "relations": "guide/FirstModel.html",
      "citation": "guide/Citation.html"
    };
    const anchor = window.location.hash.slice(1);
    if (Object.hasOwn(movedSections, anchor)) {
      const target = new URL(movedSections[anchor], window.location.href);
      target.search = window.location.search;
      target.hash = window.location.hash;
      window.location.replace(target.href);
      return;
    }
  }
  const sidebar = document.getElementById("helpSidebar");
  const scrim = document.getElementById("sidebarScrim");
  const menuButton = document.getElementById("menuButton");
  const languageButton = document.getElementById("languageButton");
  const searchInput = document.getElementById("helpSearch");
  const searchCount = document.getElementById("searchCount");
  const searchEmpty = document.getElementById("searchEmpty");
  const searchableSections = [...document.querySelectorAll("[data-searchable]")];
  const navLinks = [...document.querySelectorAll(".nav a[href^='#']")];
  const helpScript = [...document.scripts].find(script => new URL(script.src || ".", window.location.href).pathname.endsWith("/help.js"));
  const documentationRoot = helpScript
    ? new URL(".", helpScript.src)
    : new URL(".", window.location.href);

  const installDocumentationHubs = () => {
    if (!sidebar || sidebar.querySelector(".docs-hubs")) return;
    const hubs = document.createElement("nav");
    hubs.className = "docs-hubs";
    hubs.setAttribute("aria-label", "Documentation sections");
    const heading = document.createElement("p");
    heading.className = "docs-hubs-label";
    const headingChinese = document.createElement("span");
    headingChinese.className = "locale-zh";
    headingChinese.textContent = "文档分类";
    const headingEnglish = document.createElement("span");
    headingEnglish.className = "locale-en";
    headingEnglish.textContent = "Documentation";
    heading.append(headingChinese, headingEnglish);
    hubs.append(heading);

    const definitions = [
      ["G", "guide/MagneticTB.html", "指南", "Guides", "/guide/"],
      ["T", "tutorial/index.html", "教程", "Tutorials", "/tutorial/"],
      ["R", "reference/index.html", "函数参考", "Reference Pages", "/reference/"],
    ];
    for (const [code, target, zh, en, pathMarker] of definitions) {
      const link = document.createElement("a");
      link.className = "docs-hub-link";
      if (window.location.pathname.includes(pathMarker)) link.classList.add("active");
      link.href = new URL(target, documentationRoot).href;

      const marker = document.createElement("span");
      marker.className = "docs-hub-code";
      marker.textContent = code;
      const title = document.createElement("span");
      title.className = "docs-hub-title";
      const chinese = document.createElement("span");
      chinese.className = "locale-zh";
      chinese.textContent = zh;
      const english = document.createElement("span");
      english.className = "locale-en";
      english.textContent = en;
      title.append(chinese, english);
      link.append(marker, title);
      hubs.append(link);
    }
    const localNavigation = sidebar.querySelector(".nav");
    sidebar.insertBefore(hubs, localNavigation);
  };

  const katexRoot = () => {
    const sourceMarker = "/docs/user-guide/web/";
    const sourceIndex = window.location.pathname.indexOf(sourceMarker);
    if (sourceIndex >= 0) {
      const tail = window.location.pathname.slice(sourceIndex + sourceMarker.length);
      const depth = Math.max(0, tail.split("/").length - 1);
      return `${"../".repeat(depth + 3)}python/src/magnetictb/web/static/katex`;
    }
    const installedMarker = "/help/";
    const installedIndex = window.location.pathname.indexOf(installedMarker);
    if (installedIndex >= 0) {
      const tail = window.location.pathname.slice(installedIndex + installedMarker.length);
      const depth = Math.max(0, tail.split("/").length - 1);
      return `${"../".repeat(depth + 1)}static/katex`;
    }
    return "../static/katex";
  };

  const renderMath = () => {
    if (!window.katex) return;
    for (const element of document.querySelectorAll("[data-tex]")) {
      if (element.dataset.katexRendered === "true") continue;
      const source = element.textContent.trim();
      try {
        window.katex.render(source, element, {
          displayMode: true,
          output: "htmlAndMathml",
          strict: "ignore",
          throwOnError: false,
          trust: false,
        });
        element.dataset.katexRendered = "true";
        element.classList.add("katex-ready");
      } catch (_) {
        element.classList.add("katex-error");
      }
    }
  };

  const loadTheme = () => {
    const script = document.createElement("script");
    script.src = `${katexRoot()}/../theme.js?v=b608643fb0f9`;
    document.head.append(script);
  };

  const loadKatex = () => {
    if (!document.querySelector("[data-tex]")) return;
    const root = katexRoot();
    const stylesheet = document.createElement("link");
    stylesheet.rel = "stylesheet";
    stylesheet.href = `${root}/katex.min.css`;
    stylesheet.dataset.magnetictbKatex = "stylesheet";
    document.head.append(stylesheet);

    if (window.katex) {
      renderMath();
      return;
    }
    const script = document.createElement("script");
    script.src = `${root}/katex.min.js`;
    script.dataset.magnetictbKatex = "script";
    script.addEventListener("load", renderMath, { once: true });
    script.addEventListener("error", () => body.classList.add("katex-unavailable"), { once: true });
    document.head.append(script);
  };

  const resolveLocale = () => {
    const requested = new URLSearchParams(window.location.search).get("lang");
    if (requested === "zh" || requested === "en") return requested;
    try {
      const stored = window.localStorage.getItem(storageKey);
      if (stored === "zh" || stored === "en") return stored;
    } catch (_) {
      // Storage can be disabled without affecting the help page.
    }
    return "zh";
  };

  const setLocale = locale => {
    const next = locale === "en" ? "en" : "zh";
    body.dataset.locale = next;
    document.documentElement.lang = next === "en" ? "en" : "zh-CN";
    languageButton.textContent = next === "zh" ? "EN" : "中文";
    languageButton.setAttribute(
      "aria-label",
      next === "zh" ? "Switch to English" : "切换到中文",
    );
    window.dispatchEvent(new CustomEvent("magnetictb:locale", { detail: { locale: next } }));
    searchInput.placeholder = body.dataset.page === "home"
      ? (next === "zh" ? "搜索文档入口" : "Search documentation")
      : (next === "zh" ? "搜索函数、模型或错误标签" : "Search functions, models, or error tags");
    try { window.localStorage.setItem(storageKey, next); } catch (_) { /* Optional preference only. */ }
    applySearch();
  };

  const closeSidebar = () => {
    sidebar.classList.remove("open");
    scrim.classList.remove("visible");
    menuButton.setAttribute("aria-expanded", "false");
  };

  const toggleSidebar = () => {
    const open = !sidebar.classList.contains("open");
    sidebar.classList.toggle("open", open);
    scrim.classList.toggle("visible", open);
    menuButton.setAttribute("aria-expanded", String(open));
  };

  const searchableText = element => {
    const visible = element.innerText || "";
    const paired = [...element.querySelectorAll("[data-search]")]
      .map(node => node.dataset.search || "")
      .join(" ");
    return `${visible} ${paired}`.toLocaleLowerCase();
  };

  function applySearch() {
    const query = searchInput.value.trim().toLocaleLowerCase();
    let matches = 0;
    for (const section of searchableSections) {
      const matched = !query || searchableText(section).includes(query);
      section.classList.toggle("search-hidden", !matched);
      if (matched) matches += 1;
    }
    const locale = body.dataset.locale;
    searchCount.textContent = query ? String(matches) : "";
    searchCount.setAttribute(
      "aria-label",
      locale === "en" ? `${matches} matching sections` : `${matches} 个匹配章节`,
    );
    searchEmpty.classList.toggle("visible", Boolean(query) && matches === 0);
  }

  const installCopyButtons = () => {
    for (const block of document.querySelectorAll(".code-block")) {
      const pre = block.querySelector("pre");
      if (!pre || block.querySelector(".copy-button")) continue;
      const button = document.createElement("button");
      button.type = "button";
      button.className = "copy-button";
      button.textContent = "Copy";
      button.setAttribute("aria-label", "Copy code");
      button.addEventListener("click", async () => {
        try {
          await navigator.clipboard.writeText(pre.innerText);
          button.textContent = body.dataset.locale === "en" ? "Copied" : "已复制";
        } catch (_) {
          button.textContent = body.dataset.locale === "en" ? "Select text" : "请手动选择";
        }
        window.setTimeout(() => { button.textContent = "Copy"; }, 1500);
      });
      block.append(button);
    }
  };

  const observeSections = () => {
    const sections = [...document.querySelectorAll("main section[id]")];
    if (!("IntersectionObserver" in window)) return;
    const observer = new IntersectionObserver(
      entries => {
        const visible = entries
          .filter(entry => entry.isIntersecting)
          .sort((left, right) => right.intersectionRatio - left.intersectionRatio)[0];
        if (!visible) return;
        for (const link of navLinks) {
          link.classList.toggle("active", link.getAttribute("href") === `#${visible.target.id}`);
        }
      },
      { rootMargin: "-18% 0px -70% 0px", threshold: [0, .2, .5] },
    );
    for (const section of sections) observer.observe(section);
  };

  menuButton.addEventListener("click", toggleSidebar);
  scrim.addEventListener("click", closeSidebar);
  languageButton.addEventListener("click", () => {
    setLocale(body.dataset.locale === "zh" ? "en" : "zh");
  });
  searchInput.addEventListener("input", applySearch);
  for (const link of navLinks) link.addEventListener("click", closeSidebar);
  window.addEventListener("keydown", event => {
    if (event.key === "Escape") closeSidebar();
    if ((event.metaKey || event.ctrlKey) && event.key.toLocaleLowerCase() === "k") {
      event.preventDefault();
      searchInput.focus();
    }
  });

  installDocumentationHubs();
  installCopyButtons();
  observeSections();
  setLocale(resolveLocale());
  loadTheme();
  loadKatex();
})();
