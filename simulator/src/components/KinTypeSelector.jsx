const KIN_TYPES = [
  { key: 'children',      label: 'Children' },
  { key: 'siblings',      label: 'Siblings' },
  { key: 'auntsUncles',   label: 'Aunts & Uncles' },
  { key: 'cousins',       label: 'Cousins' },
  { key: 'niecesNephews', label: 'Nieces & Nephews' },
]

export default function KinTypeSelector({ selectedKin, onChange }) {
  const toggle = key => onChange({ ...selectedKin, [key]: !selectedKin[key] })

  return (
    <div className="control-group">
      <span className="control-label">Kin Types</span>
      {KIN_TYPES.map(({ key, label }) => (
        <label key={key} className="checkbox-label">
          <input
            type="checkbox"
            checked={!!selectedKin[key]}
            onChange={() => toggle(key)}
          />
          {label}
        </label>
      ))}
    </div>
  )
}
