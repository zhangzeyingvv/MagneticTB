(() => {
  "use strict";

  // Shared by the workbench and help pages; no model/API requests are made.
  const storageKey = "magnetictb.web.theme";
  const modes = ["system", "light", "dark"];
  const normalize = value => modes.includes(value) ? value : "system";
  const root = document.documentElement;
  const button = document.getElementById("themeButton");
  const system = window.matchMedia("(prefers-color-scheme: dark)");
  let preference = "system";
  try {
    preference = normalize(window.localStorage.getItem(storageKey));
  } catch (_) { /* Switching works even when saving preferences is disabled. */ }

  const updateButton = () => {
    if (!button) return;
    const english = root.lang.startsWith("en");
    const labels = english
      ? { system: "System", light: "Light", dark: "Dark" }
      : { system: "跟随系统", light: "浅色", dark: "深色" };
    const next = modes[(modes.indexOf(preference) + 1) % modes.length];
    const label = `${english ? "Theme: " : "主题："}${labels[preference]}`;
    const action = english ? `Switch to ${labels[next].toLowerCase()}` : `切换为${labels[next]}`;
    if (button.textContent !== label) button.textContent = label;
    button.title = action;
    button.setAttribute("aria-label", `${label}; ${action}`);
    button.dataset.themePreference = preference;
  };

  const apply = () => {
    root.dataset.theme = preference === "system" ? (system.matches ? "dark" : "light") : preference;
    updateButton();
  };

  if (button) {
    button.disabled = false;
    button.addEventListener("click", () => {
      preference = modes[(modes.indexOf(preference) + 1) % modes.length];
      apply();
      try { window.localStorage.setItem(storageKey, preference); } catch (_) { /* Optional storage. */ }
    });
  }
  system.addEventListener("change", () => { if (preference === "system") apply(); });
  window.addEventListener("storage", event => {
    if (event.key === storageKey || event.key === null) {
      preference = normalize(event.newValue);
      apply();
    }
  });
  window.addEventListener("magnetictb:locale", updateButton);
  // A cached page restored by Back/Forward should also use the latest setting.
  window.addEventListener("pageshow", () => {
    try { preference = normalize(window.localStorage.getItem(storageKey)); } catch (_) { /* Keep current mode. */ }
    apply();
  });
  apply();
})();
