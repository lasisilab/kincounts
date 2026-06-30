// Zero-dependency export helpers: CSV download + chart image (PNG/SVG) download.
// Everything is browser-native (Blob, URL.createObjectURL, <a download>,
// XMLSerializer, canvas) so no new packages are needed.

// Serialize an array of row objects to CSV text.
export function toCSV(rows, columns) {
  if (!rows || !rows.length) return ''
  const cols = columns ?? Object.keys(rows[0])
  const esc = v => {
    if (v == null) return ''
    const s = String(v)
    return /[",\n]/.test(s) ? `"${s.replace(/"/g, '""')}"` : s
  }
  return cols.join(',') + '\n' +
    rows.map(r => cols.map(c => esc(r[c])).join(',')).join('\n') + '\n'
}

// Trigger a browser download of arbitrary content.
export function downloadBlob(filename, content, mime = 'text/plain') {
  const blob = content instanceof Blob ? content : new Blob([content], { type: mime })
  const url = URL.createObjectURL(blob)
  const a = document.createElement('a')
  a.href = url
  a.download = filename
  document.body.appendChild(a)
  a.click()
  a.remove()
  setTimeout(() => URL.revokeObjectURL(url), 1000)
}

export function downloadCSV(filename, rows, columns) {
  downloadBlob(filename, toCSV(rows, columns), 'text/csv;charset=utf-8')
}

export function downloadText(filename, text, mime = 'text/plain;charset=utf-8') {
  downloadBlob(filename, text, mime)
}

// Grab the <svg> element a recharts ResponsiveContainer renders, from a wrapper ref.
export function svgFromRef(ref) {
  return ref?.current?.querySelector('svg') ?? null
}

function serializeSvg(svg) {
  const clone = svg.cloneNode(true)
  const rect = svg.getBoundingClientRect()
  const w = Math.ceil(rect.width) || 800
  const h = Math.ceil(rect.height) || 400
  clone.setAttribute('width', w)
  clone.setAttribute('height', h)
  clone.setAttribute('xmlns', 'http://www.w3.org/2000/svg')
  // Opaque white background so the export is not transparent.
  const bg = document.createElementNS('http://www.w3.org/2000/svg', 'rect')
  bg.setAttribute('width', '100%')
  bg.setAttribute('height', '100%')
  bg.setAttribute('fill', '#ffffff')
  clone.insertBefore(bg, clone.firstChild)
  return { xml: new XMLSerializer().serializeToString(clone), w, h }
}

export function downloadChartSvg(svg, filename) {
  if (!svg) return
  const { xml } = serializeSvg(svg)
  downloadBlob(filename, xml, 'image/svg+xml;charset=utf-8')
}

// Rasterize the SVG to PNG via an offscreen canvas (2× for crispness).
export function downloadChartPng(svg, filename, scale = 2) {
  if (!svg) return
  const { xml, w, h } = serializeSvg(svg)
  const svgUrl = 'data:image/svg+xml;charset=utf-8,' + encodeURIComponent(xml)
  const img = new Image()
  img.onload = () => {
    const canvas = document.createElement('canvas')
    canvas.width = w * scale
    canvas.height = h * scale
    const ctx = canvas.getContext('2d')
    ctx.scale(scale, scale)
    ctx.drawImage(img, 0, 0)
    canvas.toBlob(blob => { if (blob) downloadBlob(filename, blob, 'image/png') }, 'image/png')
  }
  // If rasterization fails (e.g. tainted canvas), fall back to vector SVG.
  img.onerror = () => downloadChartSvg(svg, filename.replace(/\.png$/, '.svg'))
  img.src = svgUrl
}
