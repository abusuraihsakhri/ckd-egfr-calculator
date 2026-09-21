import { calculate } from './calculator.js';

const form = document.querySelector('#calculator-form');
const results = document.querySelector('#results');
const errorBox = document.querySelector('#error-message');
const themeButton = document.querySelector('#theme-toggle');
const resetButton = document.querySelector('#reset-button');

function optionalNumber(id) {
  const raw = document.querySelector(id).value.trim();
  if (!raw) return null;
  const parsed = Number(raw);
  return Number.isFinite(parsed) ? parsed : Number.NaN;
}

function setText(id, value) {
  document.querySelector(id).textContent = value;
}

function format(value) {
  return value === null || value === undefined ? '—' : value.toFixed(1);
}

function render(data) {
  setText('#primary-egfr', format(data.primary));
  setText('#gfr-category', data.gfr);
  setText('#gfr-description', data.gfrDescription);
  setText('#albuminuria-category', data.albuminuria ?? '—');
  setText('#albuminuria-description', data.albuminuriaDescription ?? 'Add ACR for category');
  setText('#kdigo-risk', data.risk ?? '—');
  setText('#creatinine-egfr', format(data.creatinineEgfr));
  setText('#cystatin-egfr', format(data.cystatinEgfr));
  setText('#combined-egfr', format(data.combinedEgfr));
  results.hidden = false;
}

form.addEventListener('submit', (event) => {
  event.preventDefault();
  errorBox.hidden = true;
  try {
    const data = calculate({
      age: document.querySelector('#age').value,
      sex: document.querySelector('#sex').value,
      creatinine: optionalNumber('#creatinine'),
      creatinineUnit: document.querySelector('#creatinine-unit').value,
      cystatin: optionalNumber('#cystatin'),
      acr: optionalNumber('#acr'),
    });
    render(data);
  } catch (error) {
    results.hidden = true;
    errorBox.textContent = error instanceof Error ? error.message : 'Unable to calculate.';
    errorBox.hidden = false;
  }
});

resetButton.addEventListener('click', () => {
  form.reset();
  results.hidden = true;
  errorBox.hidden = true;
  document.querySelector('#age').focus();
});

const savedTheme = localStorage.getItem('theme');
if (savedTheme === 'dark' || savedTheme === 'light') document.documentElement.dataset.theme = savedTheme;

themeButton.addEventListener('click', () => {
  const dark = document.documentElement.dataset.theme === 'dark';
  const next = dark ? 'light' : 'dark';
  document.documentElement.dataset.theme = next;
  localStorage.setItem('theme', next);
  themeButton.setAttribute('aria-label', `Use ${dark ? 'dark' : 'light'} theme`);
});
