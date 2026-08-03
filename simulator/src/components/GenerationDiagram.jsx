// Compact diagram clarifying WHOSE distribution the sibling curve describes.
//
// The fertility distribution describes the women observed in a given year. The
// sibling distribution describes a DIFFERENT generation — their children — and
// is the size-biased transform of fertility. Plotting "sibling mean" against an
// observation year is otherwise easy to misread as siblings measured in that
// year, which it is not.
//
// `mothers` describes who the fertility distribution was measured on; it comes
// from the active dataset (see describeMothers in lib/datasets.js), since only
// the dataset knows whether it has a documented sample definition.
export default function GenerationDiagram({ year, mothers }) {
  return (
    <div className="gen-diagram">
      <div className="gen-diagram-row">
        <div className="gen-diagram-node mothers">
          <span className="gen-diagram-tag">Observed at {year}</span>
          <span className="gen-diagram-title">Mothers</span>
          <span className="gen-diagram-sub">{mothers}</span>
          <span className="gen-diagram-metric">completed fertility = X</span>
        </div>

        <div className="gen-diagram-arrow" aria-hidden="true">
          <span className="gen-diagram-arrow-glyph">→</span>
          <span className="gen-diagram-arrow-label">size-biasing</span>
        </div>

        <div className="gen-diagram-node children">
          <span className="gen-diagram-tag">next generation</span>
          <span className="gen-diagram-title">Their children</span>
          <span className="gen-diagram-sub">a child drawn at random from this cohort’s offspring</span>
          <span className="gen-diagram-metric">sibship size = Y</span>
        </div>
      </div>
      <p className="gen-diagram-formula">
        P(Y = k) = (k + 1) · P(X = k + 1) / E[X] — larger families contain more
        children, so they are over-represented when you sample a child rather than a mother.
      </p>
    </div>
  )
}
