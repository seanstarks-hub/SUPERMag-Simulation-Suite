# supermag.themes

Matplotlib figure styling presets. Four built-in themes for publication,
presentation, drafting, and dark environments. Supports custom theme
registration.

---

## `apply_theme()`

Apply a named theme globally by modifying `matplotlib.rcParams` in place.

```python
supermag.apply_theme(name)
```

| Name | Type | Description |
|------|------|-------------|
| `name` | str | Theme name (one of `list_themes()`) |

---

## `theme_context()`

Context manager that temporarily applies a theme, then restores defaults on exit.

```python
supermag.theme_context(name)
```

| Name | Type | Description |
|------|------|-------------|
| `name` | str | Theme name |

**Yields:** `dict` — the theme's rcParams dict.

### Example

```python
import supermag

with supermag.theme_context("publication"):
    fig, ax = plt.subplots()
    ax.plot([0, 1], [0, 1])
# rcParams restored here
```

---

## `list_themes()`

Return a sorted list of available theme names.

```python
supermag.list_themes()
```

**Returns:** `list[str]` — e.g. `['dark', 'draft', 'presentation', 'publication']`

---

## `get_theme()`

Return a copy of the rcParams dict for the named theme.

```python
supermag.get_theme(name)
```

| Name | Type | Description |
|------|------|-------------|
| `name` | str | Theme name |

**Returns:** `dict` — matplotlib rcParams overrides.

**Raises:** `KeyError` if theme name is not registered.

---

## `register_theme()`

Register a custom theme.

```python
from supermag.themes import register_theme

register_theme(name, params)
```

| Name | Type | Description |
|------|------|-------------|
| `name` | str | Theme name (overwrites existing if same name) |
| `params` | dict | Matplotlib rcParams overrides |

### Example

```python
from supermag.themes import register_theme

register_theme("minimal", {
    "figure.figsize": (4, 3),
    "font.size": 9,
    "axes.grid": False,
})
```

---

## Built-in Themes

### `publication`

APS/PRB single-column journal formatting.

| Property | Value |
|----------|-------|
| Figure size | 3.375" × 2.53" |
| Font | Serif, 8 pt |
| DPI | 300 |
| Line width | 1.0 |
| Grid | Off |
| Ticks | In-pointing, all sides |

### `presentation`

Large fonts and bold colors for talks and slides.

| Property | Value |
|----------|-------|
| Figure size | 10" × 5.625" (16:9) |
| Font | Sans-serif, 16 pt |
| DPI | 150 |
| Line width | 2.5 |
| Grid | On, α = 0.3 |
| Ticks | Out-pointing |

### `draft`

Fast iteration with low-resolution rendering.

| Property | Value |
|----------|-------|
| Figure size | 8" × 5" |
| Font | Sans-serif, 11 pt |
| DPI | 72 |
| Line width | 1.5 |
| Grid | On, dashed, α = 0.4 |
| Background | `#fafafa` |

### `dark`

Dark background for IDEs and notebook environments.

| Property | Value |
|----------|-------|
| Figure size | 8" × 5" |
| Font | Sans-serif, 12 pt |
| DPI | 100 |
| Line width | 1.8 |
| Figure background | `#1e1e1e` |
| Axes background | `#252526` |
| Text color | `#d4d4d4` |
| Grid | On, α = 0.5 |
