"""Check the shipped HTML in Chromium; save screenshots and a print-render PDF.

Run from the repository root: python docs/check_revision_browser.py
Requires Playwright and Chromium. Set CHROMIUM_EXECUTABLE for a system browser.
These tests do not constitute a formal accessibility or cross-browser audit.
"""
from __future__ import annotations
import hashlib
import json
import math
import os
import re
from pathlib import Path
from playwright.sync_api import sync_playwright


def main() -> None:
    root = Path(__file__).resolve().parents[1]
    source = root / 'stochastic_processes_revision.html'
    output = root / 'validation'
    output.mkdir(exist_ok=True)
    errors: list[str] = []
    checks: list[str] = []
    with sync_playwright() as p:
        executable = os.environ.get('CHROMIUM_EXECUTABLE')
        browser = p.chromium.launch(executable_path=executable, args=['--no-sandbox'])
        page = browser.new_page(viewport={'width': 1440, 'height': 1000}, reduced_motion='reduce')
        page.on('pageerror', lambda error: errors.append(str(error)))
        page.set_content(source.read_text(encoding="utf-8"))
        assert page.locator('details.check').count() == 20
        ids = page.locator('[id]').evaluate_all('(els) => els.map(e => e.id)')
        assert len(ids) == len(set(ids)), 'Duplicate IDs'
        for href in page.locator('a[href]').evaluate_all('(els) => els.map(e => e.getAttribute("href"))'):
            if href.startswith('#'):
                assert href[1:] in ids, href
            elif not re.match(r'^https?://', href):
                assert (root / href.split('#')[0]).is_file(), href
        checks.append('20 self-checks; unique IDs; valid internal and relative links')
        page.locator('#answers').click()
        assert page.locator('details.check[open]').count() == 20
        page.locator('#answers').click()
        assert page.locator('details.check[open]').count() == 0
        checks.append('Expand and collapse all self-checks')
        # Range controls are not text fields: use native keyboard interaction.
        slider = page.locator('#gap')
        for gap in (.05, .5, 1., 3.):
            slider.focus()
            slider.press('Home')
            for _ in range(round((gap - .05) / .05)):
                slider.press('ArrowRight')
            assert math.isclose(float(slider.input_value()), gap)
            text = page.locator('#gap-result').inner_text()
            assert f'{25 + 4 * math.exp(-gap):.4f}' in text, text
            assert f'{-2 * math.expm1(-2 * gap):.4f}' in text, text
        slider.evaluate("el => { el.value='0.5'; el.dispatchEvent(new Event('input', {bubbles:true})); }")
        checks.append('Keyboard-operated OU calculator: four gaps')
        widths = (320, 375, 768, 1024, 1440)
        for width in widths:
            page.set_viewport_size({'width': width, 'height': 1000})
            page.evaluate('window.scrollTo(0,0)')
            if page.evaluate('document.documentElement.scrollWidth > innerWidth + 1'):
                errors.append(f'Horizontal page overflow at {width}px')
            page.screenshot(path=str(output / f'viewport-{width}.png'))
        checks.append('Page layout: 320, 375, 768, 1024 and 1440 pixels')
        page.set_viewport_size({'width': 1440, 'height': 1000})
        for section in ('ar', 'brownian', 'ou', 'estimation', 'sources'):
            page.locator('#' + section).screenshot(path=str(output / f'section-{section}.png'))
        before = page.locator('details[open]').count()
        page.evaluate("window.dispatchEvent(new Event('beforeprint'))")
        assert page.locator('details[open]').count() == page.locator('details').count()
        page.evaluate("window.dispatchEvent(new Event('afterprint'))")
        assert page.locator('details[open]').count() == before
        page.pdf(path=str(output / 'revision-print.pdf'), format='A4', print_background=True)
        checks.append('Print expansion/restoration and A4 PDF render')
        context = browser.new_context(java_script_enabled=False, viewport={'width':375,'height':900}, reduced_motion='reduce')
        fallback = context.new_page()
        fallback.set_content(source.read_text(encoding="utf-8"))
        fallback.locator('details.check summary').first.click()
        assert fallback.locator('details.check[open]').count() == 1
        checks.append('No-JavaScript reading and native disclosure')
        report = {'commit':os.environ.get('GITHUB_SHA'), 'html_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),
                  'browser':browser.version, 'checks':checks, 'errors':errors}
        (output/'browser-report.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
        print(json.dumps(report,indent=2))
        browser.close()
    assert not errors, errors


if __name__ == '__main__':
    main()
