// "About this data" block — makes the source of the empirical fertility data
// explicit in the UI. Reads from the active dataset's provenance object, so it
// shows the IPUMS citation by default and the user's own description after import.
export default function DataProvenance({ dataset }) {
  const p = dataset?.provenance ?? {}
  const isUser = dataset?.id !== 'ipums'

  return (
    <details className="data-provenance" open={!isUser}>
      <summary>
        <span className="data-provenance-label">About this data</span>
        <span className="data-provenance-source">{dataset?.label ?? 'Dataset'}</span>
      </summary>
      <div className="data-provenance-body">
        {isUser && (
          <p className="data-provenance-userflag">
            Showing <strong>your imported data</strong>. It has replaced the IPUMS
            census data throughout the Fertility Fit and Simulator.
          </p>
        )}
        <dl className="data-provenance-dl">
          {p.source   && (<><dt>Source</dt><dd>{p.source}</dd></>)}
          {p.variable && (<><dt>Measure</dt><dd>{p.variable}</dd></>)}
          {p.sample   && (<><dt>Sample</dt><dd>{p.sample}</dd></>)}
          {p.cohorts  && (<><dt>Cohorts</dt><dd>{p.cohorts}</dd></>)}
        </dl>
        {p.caveat && <p className="data-provenance-caveat">{p.caveat}</p>}
        {p.methodsHref && (
          <p className="data-provenance-link">
            <a href={p.methodsHref} target="_blank" rel="noreferrer">
              {p.methodsLabel ?? 'Full methods'} ↗
            </a>
          </p>
        )}
      </div>
    </details>
  )
}
