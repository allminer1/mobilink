from click.testing import CliRunner

import linkstat
from linkstat.cli import main


def test_version_string():
    assert isinstance(linkstat.__version__, str) and linkstat.__version__


def test_cli_help_and_version():
    runner = CliRunner()
    r = runner.invoke(main, ["--help"])
    assert r.exit_code == 0 and "co-localization" in r.output
    r = runner.invoke(main, ["--version"])
    assert r.exit_code == 0 and linkstat.__version__ in r.output
