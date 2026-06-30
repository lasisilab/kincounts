import { useState } from 'react'
import { buildUserDataset, makeTemplateCSV } from '../lib/importData.js'
import { downloadText } from '../lib/exportUtils.js'

// Panel for importing your own fertility data. The app fits Poisson + ZINB to
// the raw counts you provide; the result replaces IPUMS everywhere. Two entry
// routes: upload a CSV file, or type/paste counts directly.
export default function DataImport({ baseDataset, onLoad, onClose }) {
  const [text, setText] = useState('')
  const [name, setName] = useState('My data')
  const [error, setError] = useState(null)
  const [warnings, setWarnings] = useState([])

  function handleFile(e) {
    const file = e.target.files?.[0]
    if (!file) return
    const reader = new FileReader()
    reader.onload = () => { setText(String(reader.result)); setError(null) }
    reader.readAsText(file)
  }

  function handleLoad() {
    const res = buildUserDataset(text, name.trim() || 'My data')
    if (res.error) { setError(res.error); setWarnings([]); return }
    setError(null)
    setWarnings(res.warnings ?? [])
    onLoad(res.dataset)
  }

  function handleTemplate() {
    downloadText('kincounts_data_template.csv', makeTemplateCSV(baseDataset))
  }

  return (
    <div className="import-panel">
      <div className="import-header">
        <h3>Use your own data</h3>
        <button className="import-close" onClick={onClose} aria-label="Close">×</button>
      </div>

      <p className="import-intro">
        Provide fertility as <strong>raw counts</strong> — the number of women with each
        number of children ever born, per year. The app computes the empirical
        distribution and fits Poisson and ZINB; your data then replaces the census
        data throughout. The top child count is treated as “that many or more”.
      </p>

      <ol className="import-steps">
        <li>
          <button className="import-btn-secondary" onClick={handleTemplate}>
            ↓ Download CSV template
          </button>
          <span className="import-hint">
            Pre-filled with the census cohorts as a worked example — replace the numbers with yours.
          </span>
        </li>
        <li>
          <label className="import-file-label">
            ↑ Upload a CSV
            <input type="file" accept=".csv,text/csv,text/plain" onChange={handleFile} />
          </label>
          <span className="import-hint">…or paste / edit the rows below.</span>
        </li>
      </ol>

      <label className="import-name">
        Dataset name
        <input type="text" value={name} onChange={e => setName(e.target.value)} />
      </label>

      <textarea
        className="import-textarea"
        value={text}
        onChange={e => { setText(e.target.value); setError(null) }}
        spellCheck={false}
        placeholder={'year,label,children,count\n1990,1931–1940,0,540\n1990,1931–1940,1,820\n1990,1931–1940,2,1980\n...'}
        rows={10}
      />

      {error && <p className="import-error">{error}</p>}
      {warnings.length > 0 && (
        <ul className="import-warnings">
          {warnings.map((w, i) => <li key={i}>{w}</li>)}
        </ul>
      )}

      <div className="import-actions">
        <button className="import-btn-primary" onClick={handleLoad} disabled={!text.trim()}>
          Load data → replace census data
        </button>
        <button className="import-btn-secondary" onClick={onClose}>Cancel</button>
      </div>
    </div>
  )
}
