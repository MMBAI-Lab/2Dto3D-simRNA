-- DansLab report styling for the pandoc DOCX writer.
--
-- The Markdown reports stay plain so they read well on their own; this
-- filter is what adapts them to the page:
--   * "**DANSLAB · MMBAI**" becomes the red mono kicker, the rest of that
--     line muted mono metadata
--   * the language-switch line becomes muted mono metadata
--   * dot-brackets get their own smaller mono style, and any table holding
--     one gets column widths wide enough to keep it on a single line
--   * structure figures are capped so a figure and its caption share a page
-- The character styles are defined in the generated reference DOCX
-- (see scripts/build_report.py).

local stringify = pandoc.utils.stringify

-- Page geometry and type metrics, in points. A4 minus 16 mm side margins.
local TEXT_WIDTH_PT = 504
local CHAR_DOTBRACKET = 4.0   -- IBM Plex Mono at 6.5 pt
local CHAR_MONO = 4.9         -- IBM Plex Mono at 8 pt
local CHAR_SANS = 5.0         -- IBM Plex Sans at 8.5 pt, caps-heavy names
local CELL_PADDING_PT = 10
local PROSE_CAP_CHARS = 26    -- prose wraps, so cap its claim on the width
local DOTBRACKET_MIN = 40     -- shorter runs fit anywhere; leave them alone
local FIGURE_WIDTH = "4.6in"

local function ulen(s)
  return utf8 and utf8.len(s) or #s
end

local function span(inlines, name)
  return pandoc.Span(inlines, { ["custom-style"] = name })
end

local function is_dotbracket(text)
  return ulen(text) >= DOTBRACKET_MIN and text:match("^[%(%)%.%[%]{}<>|_,%-]+$") ~= nil
end

-- Report header boilerplate -------------------------------------------------

function Para(el)
  local first = el.content[1]
  if first and first.t == "Strong" and stringify(first):match("^DANSLAB") then
    local rest = {}
    for i = 2, #el.content do
      rest[#rest + 1] = el.content[i]
    end
    return pandoc.Para({ span(first.content, "Kicker"), span(rest, "Meta") })
  end

  local is_lang_switch = false
  pandoc.walk_block(el, {
    Link = function(link)
      if link.target:match("^REPORT_E[NS]%.md$") or link.target:match("^README_E[NS]%.md$") then
        is_lang_switch = true
      end
    end,
  })
  if is_lang_switch then
    return pandoc.Para({ span(el.content, "Meta") })
  end
  return nil
end

-- Dot-brackets --------------------------------------------------------------

function Code(el)
  if is_dotbracket(el.text) then
    return span({ pandoc.Str(el.text) }, "DotBracket")
  end
  return nil
end

-- Tables --------------------------------------------------------------------

-- Width a cell needs, in points. Unbreakable runs (dot-brackets, long code)
-- claim their full length; prose is capped because it wraps without harm.
-- Measured on the flattened text: Code inlines have already been rewritten
-- into DotBracket spans by the time block filters run.
-- Returns: width the cell wants, width below which it starts breaking words.
local function cell_demand(cell)
  local text = stringify(pandoc.Div(cell))
  local rigid, longest_word = 0, 0
  for token in text:gmatch("%S+") do
    local len = ulen(token)
    if len >= DOTBRACKET_MIN then
      local per_char = is_dotbracket(token) and CHAR_DOTBRACKET or CHAR_MONO
      rigid = math.max(rigid, len * per_char)
    else
      longest_word = math.max(longest_word, len * CHAR_SANS)
    end
  end
  local want = math.max(rigid, longest_word,
                        math.min(ulen(text), PROSE_CAP_CHARS) * CHAR_SANS)
  local floor = math.max(rigid, longest_word)
  return want + CELL_PADDING_PT, floor + CELL_PADDING_PT, rigid > 0
end

-- Column widths as fractions of the text width. Every column keeps at least
-- the width of its longest unbreakable run — a dot-bracket or a tool name —
-- and whatever room is left over goes to the columns that asked for more.
-- Only tables containing a long dot-bracket are sized; the rest are left to
-- Word's auto-fit, which handles prose better than any estimate here.
local function column_widths(header, rows, ncols)
  local wants, floors, any_rigid = {}, {}, false
  for i = 1, ncols do
    wants[i], floors[i] = 0, 0
  end
  local function scan(row)
    for i, cell in ipairs(row) do
      if i <= ncols then
        local want, floor, rigid = cell_demand(cell)
        wants[i] = math.max(wants[i], want)
        floors[i] = math.max(floors[i], floor)
        any_rigid = any_rigid or rigid
      end
    end
  end
  if header then scan(header) end
  for _, row in ipairs(rows) do scan(row) end
  if not any_rigid then
    return nil
  end

  local floor_total, slack_total = 0, 0
  for i = 1, ncols do
    floor_total = floor_total + floors[i]
    slack_total = slack_total + math.max(0, wants[i] - floors[i])
  end

  local widths, total = {}, 0
  local room = TEXT_WIDTH_PT - floor_total
  for i = 1, ncols do
    local w = floors[i]
    if room > 0 then
      if slack_total > 0 then
        w = w + room * math.max(0, wants[i] - floors[i]) / slack_total
      else
        w = w + room / ncols
      end
    end
    widths[i] = w
    total = total + w
  end
  for i = 1, ncols do
    widths[i] = widths[i] / total  -- the table itself is written at 100% width
  end
  return widths
end

function Table(el)
  -- pandoc <= 2.9 simple-table AST
  if el.widths and el.headers and el.rows then
    local widths = column_widths(el.headers, el.rows, #el.widths)
    if widths then
      el.widths = widths
      return el
    end
    return nil
  end

  -- pandoc >= 2.10 AST: colspecs plus head/bodies
  if el.colspecs then
    local rows = {}
    local function collect(row_list)
      for _, row in ipairs(row_list) do
        local cells = {}
        for _, cell in ipairs(row.cells) do
          cells[#cells + 1] = cell.contents
        end
        rows[#rows + 1] = cells
      end
    end
    if el.head then collect(el.head.rows) end
    for _, body in ipairs(el.bodies or {}) do collect(body.body) end
    local widths = column_widths(nil, rows, #el.colspecs)
    if widths then
      for i, spec in ipairs(el.colspecs) do
        el.colspecs[i] = { spec[1], widths[i] }
      end
      return el
    end
  end
  return nil
end

-- Figures -------------------------------------------------------------------

function Image(el)
  if el.src:match("danslab%-logo") then
    return nil  -- the logo is sized by its own dpi metadata
  end
  if not el.attributes.width and not el.attributes.height then
    el.attributes.width = FIGURE_WIDTH
    return el
  end
  return nil
end
