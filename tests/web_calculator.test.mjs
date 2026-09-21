import test from 'node:test';
import assert from 'node:assert/strict';
import {
  egfrCreatinine2021, egfrCystatin2012, egfrCombined2021,
  gfrCategory, albuminuriaCategory, calculate
} from '../calculator.js';

test('2021 creatinine reference case', () => {
  assert.ok(Math.abs(egfrCreatinine2021(60, 'M', 1.2) - 69.2) < 0.5);
});

test('2012 cystatin C equation uses 0.996 age factor', () => {
  const expected = 133 * (0.996 ** 50) * 0.932;
  assert.ok(Math.abs(egfrCystatin2012(50, 'F', 0.8) - expected) < 1e-10);
});

test('combined equation matches formula', () => {
  const expected = 135 * ((1.3 / 0.7) ** -0.544) * ((1.1 / 0.8) ** -0.778) * (0.9961 ** 65) * 0.963;
  assert.ok(Math.abs(egfrCombined2021(65, 'F', 1.3, 1.1) - expected) < 1e-10);
});

test('GFR boundary gaps are closed', () => {
  assert.equal(gfrCategory(89.9999995)[0], 'G2');
  assert.equal(gfrCategory(59.9999995)[0], 'G3a');
  assert.equal(gfrCategory(14.9999995)[0], 'G5');
});

test('albuminuria thresholds follow KDIGO categories', () => {
  assert.equal(albuminuriaCategory(29.999)[0], 'A1');
  assert.equal(albuminuriaCategory(30)[0], 'A2');
  assert.equal(albuminuriaCategory(300)[0], 'A3');
});

test('browser calculation prefers combined equation when both markers are present', () => {
  const result = calculate({ age: 65, sex: 'F', creatinine: 1.3, creatinineUnit: 'mg/dL', cystatin: 1.1, acr: 150 });
  assert.equal(result.primary, result.combinedEgfr);
  assert.equal(result.albuminuria, 'A2');
});
