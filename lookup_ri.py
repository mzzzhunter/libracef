"""Look up mean temperature-programmed WAX and 5MS retention indices.

Usage: python lookup_ri.py 67-56-1
Missing values are empty strings; network and page-format errors are explicit.
"""

import argparse
import json
import math
import re
import sys
from statistics import fmean
from urllib.error import HTTPError, URLError
from urllib.parse import urlencode
from urllib.request import Request, urlopen

from bs4 import BeautifulSoup


class LookupRIError(RuntimeError):
    """NIST could not be reached or its response could not be interpreted."""


def fetch_ri_html(cas_number: str) -> bytes:
    """Validate a hyphenated CAS number and fetch its WebBook GC page."""
    if not isinstance(cas_number, str) or not re.fullmatch(r'[0-9]{2,7}-[0-9]{2}-[0-9]', cas_number.strip()):
        raise ValueError('CAS number must have the form 67-56-1.')
    query = urlencode({'ID': cas_number.strip(), 'Units': 'SI', 'Mask': '2000'})
    request = Request('https://webbook.nist.gov/cgi/cbook.cgi?' + query,
                      headers={'User-Agent': 'LibraCEF/0.4 (RI lookup)'})
    try:
        with urlopen(request, timeout=30) as response:
            return response.read()
    except HTTPError as error:
        raise LookupRIError(f'NIST WebBook returned HTTP {error.code}: {error.reason}') from error
    except (URLError, OSError) as error:
        raise LookupRIError(f'Could not read NIST WebBook: {error}') from error


def _text(element) -> str:
    return ' '.join(element.stripped_strings).strip()


def _phase_category(phase: str) -> str | None:
    phase = re.sub(r'\s*-\s*', '-', phase.strip().upper())
    if re.fullmatch(r'[A-Z0-9]+(?:-[A-Z0-9]+)*-WAX', phase):
        return 'WAX'
    if phase in ('HP-5MS', 'DB-5MS'):
        return '5MS'
    return None


def parse_ri_html(html: str | bytes) -> dict[str, float | str]:
    """Average matching numeric RI rows inside the gas chromatography section.

    Every reported row has equal weight. Ramp and custom programmes are
    included, even when a custom programme's comment omits its details.
    """
    soup = BeautifulSoup(html, 'html.parser')
    values = {'WAX': [], '5MS': []}
    section = next((h for h in soup.find_all('h2')
                    if _text(h).casefold() == 'gas chromatography'), None)
    if section is None:
        text = soup.get_text(' ', strip=True).casefold()
        if 'cas registry number:' in text or 'species not found' in text or 'not found' in [
                _text(h).casefold() for h in soup.find_all('h1')]:
            return {'WAX': '', '5MS': ''}
        raise LookupRIError('Unrecognised NIST WebBook page: no compound or chromatography section.')

    heading = ''
    for element in section.find_all_next(['h2', 'h3', 'table']):
        if element.name == 'h2':
            break
        if element.name == 'h3':
            heading = _text(element).casefold()
            continue
        if 'isothermal' in heading:
            continue
        if not any(term in heading for term in ('temperature ramp', 'temperature program', 'temperature-program')):
            raise LookupRIError(f'Unrecognised chromatography table conditions: {heading!r}')
        rows = element.find_all('tr')
        if not rows:
            raise LookupRIError('Chromatography table has no column headers.')
        headers = [_text(cell).casefold() for cell in rows[0].find_all(['th', 'td'])]
        if 'active phase' not in headers or 'i' not in headers:
            raise LookupRIError('Chromatography table is missing Active phase or I columns.')
        phase_index, ri_index = headers.index('active phase'), headers.index('i')
        for row in rows[1:]:
            cells = row.find_all(['td', 'th'], recursive=False)
            if not cells:
                continue
            if len(cells) != len(headers):
                raise LookupRIError('Chromatography row does not match its column headers.')
            category = _phase_category(_text(cells[phase_index]))
            if category is None:
                continue
            try:
                ri = float(_text(cells[ri_index]))
            except ValueError:
                continue
            if math.isfinite(ri):
                values[category].append(ri)
    return {category: fmean(entries) if entries else '' for category, entries in values.items()}


def lookup_ri(cas_number: str) -> dict[str, float | str]:
    """Return {'WAX': mean or '', '5MS': mean or ''} for a CAS number."""
    return parse_ri_html(fetch_ri_html(cas_number))


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('cas_number', help='Hyphenated CAS number, e.g. 67-56-1')
    args = parser.parse_args(argv)
    try:
        result = lookup_ri(args.cas_number)
    except (ValueError, LookupRIError) as error:
        print(f'RI lookup failed: {error}', file=sys.stderr)
        return 1
    print(json.dumps(result, allow_nan=False))
    return 0


if __name__ == '__main__':
    sys.exit(main())
