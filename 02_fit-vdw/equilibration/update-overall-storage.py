"""
Script 
"""
import logging
import tempfile

import click
from eveq.storage.storage import LocalStoredEquilibrationData

logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")

@click.command()
@click.option(
    "--lfs-root-directory",
    "-lfs",
    "lfs_root_directory",
    type=click.Path(exists=True, dir_okay=True, readable=True),
    default="stored_data",
    help=(
        "Path to the root directory of an existing LocalFileStorage. "
        "The default path here is from the equilibration of properties not already found in ../equilibration_stored_data."
    )
)
@click.option(
    "--lsed-root-directory",
    "-lsed",
    "lsed_root_directory",
    type=click.Path(file_okay=False, writable=True),
    default="../equilibration_stored_data",
    help=(
        "Path to the ongoing directory for the LocalStoredEquilibrationData. "
        "This was created from equilibrating ash-sage-rc1."
    )
)
def main(
    lfs_root_directory: str = "stored_data",
    lsed_root_directory: str = "../equilibration_stored_data",
):
    temporary_directory = tempfile.mkdtemp()
    logger.info(f"Temporary directory created at {temporary_directory}")

    storage = LocalStoredEquilibrationData.from_localfilestorage(
        lfs_root_directory=lfs_root_directory,
        new_root_directory=str(temporary_directory)
    )
    logger.info(f"Loaded {len(storage._cached_retrieved_objects)} new objects from LocalFileStorage")

    ongoing_storage = LocalStoredEquilibrationData(lsed_root_directory)
    logger.info(f"Loaded {len(ongoing_storage._cached_retrieved_objects)} objects in existing storage")

    # Combine the two storages
    ongoing_storage.update(storage)
    logger.info(f"Combined storage now has {len(ongoing_storage._cached_retrieved_objects)} objects")


if __name__ == "__main__":
    main()
