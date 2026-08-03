import { useState } from 'react'
import { buildUserDataset } from '../lib/importData.js'
import {
  makeTemplateCSV,
  poissonExampleCSV,
  overdispersedExampleCSV,
  POISSON_EXAMPLE,
  OVERDISPERSED_EXAMPLE,
} from '../lib/exampleData.js'
import { downloadText } from '../lib/exportUtils.js'

// Panel for importing your own fertility data. The app fits Poisson + ZINB to
// the raw counts you provide; the result replaces IPUMS everywhere. Two entry
// routes: upload a CSV file, or type/paste counts directly.
export default function DataImport({ onLoad, onClose }) {
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
    downloadText('kincounts_data_template.csv', makeTemplateCSV())
  }

  // Drop a synthetic example into the textarea rather than loading it straight
  // in, so the counts can be inspected (and edited) before they are fitted.
  function loadExample(csv, exampleName) {
    setText(csv)
    setName(exampleName)
    setError(null)
    setWarnings([])
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
        distribution and fits Poisson and ZINB models to your data.
        The top child count is treated as “that many or more”.
      </p>

      <ol className="import-steps">
        <li>
          <button className="import-btn-secondary" onClick={handleTemplate}>
            ↓ Download CSV template
          </button>
          <span className="import-hint">
            Two placeholder cohorts showing the required columns — replace the numbers with yours.
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

      <div className="import-examples">
        <span className="import-examples-label">Or load a synthetic example</span>
        <button
          className="import-btn-secondary"
          onClick={() => loadExample(poissonExampleCSV(), POISSON_EXAMPLE.name)}
        >
          Poisson (λ = 1.8, 2.2)
        </button>
        <button
          className="import-btn-secondary"
          onClick={() => loadExample(overdispersedExampleCSV(), OVERDISPERSED_EXAMPLE.name)}
        >
          Overdispersed (θ = 2.5, 3.0)
        </button>
        <span className="import-hint">
          Generated from the named distribution — illustrative, not real data. Because the
          truth is known, these check the fit: the Poisson example should select
          <strong> Poisson</strong> (ZINB can only match it, and pays for two extra
          parameters), the overdispersed one should select <strong>ZINB</strong> (Poisson
          cannot represent variance above the mean).
        </span>
      </div>

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
          Load data
        </button>
        <button className="import-btn-secondary" onClick={onClose}>Cancel</button>
      </div>
    </div>
  )
}
