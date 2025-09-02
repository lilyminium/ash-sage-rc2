import pathlib
import json
import click
from loguru import logger

@click.command()
@click.option(
    "--name-and-scope-key-file",
    multiple=True,
    type=(str, str),
    help="Name and scope key file pairs"
)
def main(
    name_and_scope_key_files: list[tuple[str, str]] = [],
    output_file: str = "scope-key-config.json"
):
    scope_keys: dict[str, str] = {}
    for name, scope_key_path in name_and_scope_key_files:
        if not pathlib.Path(scope_key_path).exists():
            raise FileNotFoundError(f"Scope key file for {name} not found at {scope_key_path}")
        scope_key = pathlib.Path(scope_key_path).read_text().strip()
        scope_keys[name] = scope_key

    with open(output_file, 'w') as f:
        json.dump({"scope_keys": scope_keys}, f, indent=4)

    logger.info(f"Scope key configuration written to {output_file}")

if __name__ == "__main__":
    main()
