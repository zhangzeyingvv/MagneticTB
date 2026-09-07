(() => {
  "use strict";

  const STORAGE_KEY = "magnetictb.web.locale";
  const SUPPORTED = new Set(["zh", "en"]);
  const ATTRIBUTES = ["aria-label", "placeholder", "title"];
  const textSources = new WeakMap();
  const attributeSources = new WeakMap();
  let locale = "zh";
  let observer = null;

  const translations = () => window.MagneticTBTranslations?.en || {};
  const replacements = () => Object.entries(translations()).sort((left, right) => right[0].length - left[0].length);
  const hasHan = value => /[\u3400-\u9fff]/u.test(value);

  const translate = value => {
    if (locale !== "en" || typeof value !== "string" || !hasHan(value)) return value;
    let translated = value;
    for (const [source, target] of replacements()) translated = translated.split(source).join(target);
    return translated;
  };

  const localizeText = node => {
    if (!textSources.has(node)) textSources.set(node, node.nodeValue || "");
    const localized = locale === "en" ? translate(textSources.get(node)) : textSources.get(node);
    if (node.nodeValue !== localized) node.nodeValue = localized;
  };

  const localizeAttributes = element => {
    let sources = attributeSources.get(element);
    if (!sources) {
      sources = {};
      for (const name of ATTRIBUTES) if (element.hasAttribute(name)) sources[name] = element.getAttribute(name);
      attributeSources.set(element, sources);
    }
    for (const [name, source] of Object.entries(sources)) {
      const localized = locale === "en" ? translate(source) : source;
      if (element.getAttribute(name) !== localized) element.setAttribute(name, localized);
    }
  };

  const localizeTree = root => {
    if (!root) return;
    if (root.nodeType === Node.TEXT_NODE) localizeText(root);
    if (root.nodeType === Node.ELEMENT_NODE) localizeAttributes(root);
    const walker = document.createTreeWalker(root, NodeFilter.SHOW_ELEMENT | NodeFilter.SHOW_TEXT);
    let node;
    while ((node = walker.nextNode())) {
      if (node.nodeType === Node.TEXT_NODE) localizeText(node);
      else localizeAttributes(node);
    }
  };

  const updateLanguageButton = () => {
    const button = document.getElementById("languageButton");
    if (!button) return;
    const label = locale === "zh" ? "EN" : "中文";
    const accessibleLabel = locale === "zh" ? "Switch to English" : "切换到中文";
    if (button.textContent !== label) button.textContent = label;
    if (button.getAttribute("aria-label") !== accessibleLabel) button.setAttribute("aria-label", accessibleLabel);
    if (button.getAttribute("title") !== accessibleLabel) button.setAttribute("title", accessibleLabel);
  };

  const apply = () => {
    document.documentElement.lang = locale === "en" ? "en" : "zh-CN";
    localizeTree(document.body);
    updateLanguageButton();
    window.dispatchEvent(new CustomEvent("magnetictb:locale", { detail: { locale } }));
  };

  const setLocale = requested => {
    locale = SUPPORTED.has(requested) ? requested : "zh";
    try { window.localStorage.setItem(STORAGE_KEY, locale); } catch (_) { /* Storage may be disabled. */ }
    apply();
    return locale;
  };

  const resolveInitialLocale = () => {
    const query = new URLSearchParams(window.location.search).get("lang");
    if (SUPPORTED.has(query)) return query;
    try {
      const stored = window.localStorage.getItem(STORAGE_KEY);
      if (SUPPORTED.has(stored)) return stored;
    } catch (_) { /* Storage may be disabled. */ }
    return "zh";
  };

  const start = () => {
    locale = resolveInitialLocale();
    apply();
    if (!observer) {
      observer = new MutationObserver(records => {
        for (const record of records) {
          if (record.type === "characterData") localizeText(record.target);
          for (const node of record.addedNodes) localizeTree(node);
        }
        updateLanguageButton();
      });
      observer.observe(document.body, { childList: true, subtree: true, characterData: true });
    }
    return locale;
  };

  window.MagneticTBI18n = Object.freeze({
    get locale() { return locale; },
    start,
    setLocale,
    toggle: () => setLocale(locale === "zh" ? "en" : "zh"),
    translate,
    localizeTree,
  });
})();
