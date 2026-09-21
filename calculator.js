export const MGDL_TO_UMOLL = 88.4;

function finiteNumber(value, label) {
  const number = Number(value);
  if (!Number.isFinite(number)) throw new Error(`${label} must be a finite number.`);
  return number;
}

export function creatinineToMgDl(value, unit) {
  const number = finiteNumber(value, "Creatinine");
  if (unit === "mg/dL") return number;
  if (unit === "µmol/L") return number / MGDL_TO_UMOLL;
  throw new Error("Unsupported creatinine unit.");
}

export function validateAdultInputs({ age, sex, creatinineMgDl = null, cystatin = null, acr = null }) {
  const parsedAge = finiteNumber(age, "Age");
  if (parsedAge < 18) throw new Error("CKD-EPI adult equations require age 18 years or older.");
  if (!['F', 'M'].includes(sex)) throw new Error("Select female or male sex for the CKD-EPI equation.");
  if (creatinineMgDl !== null && (!Number.isFinite(creatinineMgDl) || creatinineMgDl <= 0)) {
    throw new Error("Creatinine must be greater than 0.");
  }
  if (cystatin !== null && (!Number.isFinite(cystatin) || cystatin <= 0)) {
    throw new Error("Cystatin C must be greater than 0.");
  }
  if (acr !== null && (!Number.isFinite(acr) || acr < 0)) throw new Error("ACR cannot be negative.");
  return parsedAge;
}

export function egfrCreatinine2021(age, sex, scr) {
  validateAdultInputs({ age, sex, creatinineMgDl: scr });
  const [kappa, alpha, sexFactor] = sex === 'F' ? [0.7, -0.241, 1.012] : [0.9, -0.302, 1];
  const ratio = scr / kappa;
  return 142 * Math.min(ratio, 1) ** alpha * Math.max(ratio, 1) ** -1.2 * 0.9938 ** age * sexFactor;
}

export function egfrCystatin2012(age, sex, scys) {
  validateAdultInputs({ age, sex, cystatin: scys });
  const ratio = scys / 0.8;
  const sexFactor = sex === 'F' ? 0.932 : 1;
  return 133 * Math.min(ratio, 1) ** -0.499 * Math.max(ratio, 1) ** -1.328 * 0.996 ** age * sexFactor;
}

export function egfrCombined2021(age, sex, scr, scys) {
  validateAdultInputs({ age, sex, creatinineMgDl: scr, cystatin: scys });
  const [kappa, alpha, sexFactor] = sex === 'F' ? [0.7, -0.219, 0.963] : [0.9, -0.144, 1];
  const scrRatio = scr / kappa;
  const scysRatio = scys / 0.8;
  return 135 * Math.min(scrRatio, 1) ** alpha * Math.max(scrRatio, 1) ** -0.544
    * Math.min(scysRatio, 1) ** -0.323 * Math.max(scysRatio, 1) ** -0.778
    * 0.9961 ** age * sexFactor;
}

export function gfrCategory(value) {
  const egfr = finiteNumber(value, "eGFR");
  if (egfr < 0) throw new Error("eGFR cannot be negative.");
  if (egfr >= 90) return ["G1", "Normal or high"];
  if (egfr >= 60) return ["G2", "Mildly decreased"];
  if (egfr >= 45) return ["G3a", "Mildly to moderately decreased"];
  if (egfr >= 30) return ["G3b", "Moderately to severely decreased"];
  if (egfr >= 15) return ["G4", "Severely decreased"];
  return ["G5", "Kidney failure"];
}

export function albuminuriaCategory(value) {
  const acr = finiteNumber(value, "ACR");
  if (acr < 0) throw new Error("ACR cannot be negative.");
  if (acr < 30) return ["A1", "Normal to mildly increased"];
  if (acr < 300) return ["A2", "Moderately increased"];
  return ["A3", "Severely increased"];
}

const riskGrid = {
  G1: { A1: "Low risk", A2: "Moderately increased risk", A3: "High risk" },
  G2: { A1: "Low risk", A2: "Moderately increased risk", A3: "High risk" },
  G3a: { A1: "Moderately increased risk", A2: "High risk", A3: "Very high risk" },
  G3b: { A1: "High risk", A2: "Very high risk", A3: "Very high risk" },
  G4: { A1: "Very high risk", A2: "Very high risk", A3: "Very high risk" },
  G5: { A1: "Very high risk", A2: "Very high risk", A3: "Very high risk" },
};

export function kdigoRisk(g, a) {
  const value = riskGrid[g]?.[a];
  if (!value) throw new Error("Unknown KDIGO category pair.");
  return value;
}

export function calculate({ age, sex, creatinine, creatinineUnit, cystatin, acr }) {
  const parsedAge = finiteNumber(age, "Age");
  const scr = creatinine === null ? null : creatinineToMgDl(creatinine, creatinineUnit);
  const scys = cystatin === null ? null : finiteNumber(cystatin, "Cystatin C");
  const parsedAcr = acr === null ? null : finiteNumber(acr, "ACR");
  if (scr === null && scys === null) throw new Error("Enter creatinine, cystatin C, or both.");
  validateAdultInputs({ age: parsedAge, sex, creatinineMgDl: scr, cystatin: scys, acr: parsedAcr });

  const creatinineEgfr = scr === null ? null : egfrCreatinine2021(parsedAge, sex, scr);
  const cystatinEgfr = scys === null ? null : egfrCystatin2012(parsedAge, sex, scys);
  const combinedEgfr = scr !== null && scys !== null ? egfrCombined2021(parsedAge, sex, scr, scys) : null;
  const primary = combinedEgfr ?? cystatinEgfr ?? creatinineEgfr;
  const [gfr, gfrDescription] = gfrCategory(primary);
  const albumin = parsedAcr === null ? null : albuminuriaCategory(parsedAcr);
  return {
    creatinineEgfr, cystatinEgfr, combinedEgfr, primary,
    gfr, gfrDescription,
    albuminuria: albumin?.[0] ?? null,
    albuminuriaDescription: albumin?.[1] ?? null,
    risk: albumin ? kdigoRisk(gfr, albumin[0]) : null,
  };
}
