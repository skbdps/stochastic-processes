"""Smoke-check the delivered standalone HTML, without external page resources.

This is a Chromium functional/layout check, not a cross-browser or full
accessibility certification. The checked file must match the edition hash.
"""
from __future__ import annotations
import hashlib
import json
import os
from pathlib import Path
from playwright.sync_api import sync_playwright
from build_notes import EXPECTED_HTML

HERE = Path(__file__).resolve().parent
SOURCE = HERE.parents[1] / 'stochastic_processes_expanded_notes.html'
OUT = HERE / 'validation'


def main() -> None:
    OUT.mkdir(exist_ok=True)
    digest = hashlib.sha256(SOURCE.read_bytes()).hexdigest()
    assert digest == EXPECTED_HTML, 'Wrong teaching edition'
    errors: list[str] = []
    requests: list[str] = []
    widths = (320, 375, 768, 1024, 1440)
    with sync_playwright() as pw:
        executable = os.getenv('CHROMIUM_EXECUTABLE')
        browser = pw.chromium.launch(**({'executable_path': executable} if executable else {}))
        page = browser.new_page(viewport={'width': 1440, 'height': 1000}, reduced_motion='reduce')
        page.on('pageerror', lambda error: errors.append(str(error)))
        page.on('request', lambda request: requests.append(request.url))
        page.goto(SOURCE.as_uri())
        assert page.locator('section.chapter').count() == 13
        assert page.locator('section.chapter h3').count() == 149
        assert page.locator('math').count() == 2147
        assert page.locator('math[display="block"]').count() == 424
        assert page.locator('math merror').count() == 0
        ids = page.locator('[id]').evaluate_all('(els) => els.map(e => e.id)')
        assert len(ids) == len(set(ids)), 'Duplicate IDs'
        for href in page.locator('a[href^="#"]').evaluate_all('(els) => els.map(e => e.getAttribute("href"))'):
            assert href[1:] in ids, href
        assert page.locator('script[src], link[rel="stylesheet"], img[src]').count() == 0
        for width in widths:
            page.set_viewport_size({'width': width, 'height': 1000})
            page.evaluate('window.scrollTo(0, 0)')
            assert page.evaluate('document.documentElement.scrollWidth <= innerWidth + 1'), width
            page.screenshot(path=str(OUT / f'viewport-{width}.png'))
        page.set_viewport_size({'width': 375, 'height': 900})
        page.locator('#open-contents').click()
        assert page.locator('#contents-disclosure').evaluate('(e) => e.open')
        page.locator('.sidebar a[href="#likelihood"]').click()
        assert page.url.endswith('#likelihood')
        page.locator('#copy-code').click()
        page.wait_for_function("['Code copied', 'Select the code below to copy'].includes(document.getElementById('copy-code').textContent)")
        page.emulate_media(media='print')
        assert not page.locator('.sidebar').is_visible()
        assert page.locator('section.chapter').count() == 13
        page.emulate_media(media='screen')
        fallback = browser.new_context(java_script_enabled=False, viewport={'width': 375, 'height': 900})
        nojs = fallback.new_page()
        nojs.goto(SOURCE.as_uri())
        assert nojs.locator('section.chapter h3').count() == 149
        assert nojs.locator('math').count() == 2147
        assert not nojs.locator('#print-book').is_visible()
        nojs.locator('#contents-disclosure summary').click()
        assert not nojs.locator('#contents-disclosure').evaluate('(e) => e.open')
        assert not errors, errors
        assert all(url.startswith('file:') for url in requests), requests
        report = {'html_sha256': digest, 'browser': browser.version, 'chapters': 13,
                  'core_subsections': 149, 'math_expressions': 2147, 'display_equations': 424,
                  'viewport_widths': widths, 'javascript_errors': errors,
                  'network_page_resources': 0, 'no_javascript_fallback': 'passed',
                  'scope': 'Standalone Chromium smoke checks; not exhaustive browser/accessibility certification.'}
        (OUT / 'browser_report.json').write_text(json.dumps(report, indent=2), encoding='utf-8')
        print(json.dumps(report, indent=2))
        browser.close()


if __name__ == '__main__':
    main()
