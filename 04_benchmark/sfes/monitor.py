import click
import os
import pathlib

from loguru import logger
from alchemiscale import AlchemiscaleClient, Scope, ScopedKey


@click.command
@click.option(
    '--scope_key',
    type=click.Path(dir_okay=False, file_okay=True, path_type=pathlib.Path),
    required=True,
    default="scoped-key.dat",
    help="Path to a serialized ScopedKey",
)
@click.option(
    '--restart',
    is_flag=True,
    help="Restart any failures",
)
@click.option(
    '--user_id',
    type=str,
    required=False,
    default=None,
)
@click.option(
    '--user_key',
    type=str,
    required=False,
    default=None,
)
@click.option(
    '--priority',
    type=float,
    required=False,
    default=None,
    help="Set the priority of the network",
)
@click.option(
    '--cancel',
    is_flag=True,
    help="Cancel all tasks in the network",
)
def run(
    scope_key: pathlib.Path,
    restart: bool,
    user_id: str,
    user_key: str,
    priority: int | None,
    cancel: bool,
):
    """
    Monitor a network of transformation tasks and
    optionally restart failed tasks.

    Parameters
    ----------
    scope_key : pathlib.Path
      A path to a serialized ScopeKey
    restart : bool
      Whether or not to attempt to restart failed tasks.
    user_id : Optional[str]
      A string for a user ID, if undefined will
      fetch from the environment variable ALCHEMISCALE_ID.
    user_key : Optional[str]
      A string for the user key, if underfined will
      fetch from the environment variable ALCHEMISCALE_KEY.
    """
    # Get the alchemiscale bits
    if user_id is None:
        user_id = os.environ['ALCHEMISCALE_ID']
    if user_key is None:
        user_key = os.environ['ALCHEMISCALE_KEY']
    asc = AlchemiscaleClient(
        'https://api.alchemiscale.org',
        user_id,
        user_key
    )

    with open(scope_key, 'r') as f:
        network_sk = f.read()

    asc.get_network_status(network_sk)

    if priority is not None:
        asc.set_network_weight(network_sk, priority)
        logger.info("Network priority set to", asc.get_network_weight(network_sk))

    if restart:
        err_tasks = asc.get_network_tasks(network_sk, status="error")
        logger.info(f"Number of errored tasks found: {len(err_tasks)}")
        if len(err_tasks) > 0:
            logger.info("Will attempt to restart tasks")
            asc.set_tasks_status(err_tasks, 'waiting')

            logger.info("Printing new network status")
            asc.get_network_status(network_sk)
        else:
            logger.info("No errored tasks were found, no further action")

    if cancel:
        task_sks = asc.get_network_tasks(network_sk)
        logger.info("Cancelling all tasks")
        asc.cancel_tasks(task_sks, network_sk)
        logger.info("Setting all tasks to deleted")
        asc.set_tasks_status(task_sks, 'deleted')
        logger.info("Printing new network status")
        asc.get_network_status(network_sk)




if __name__ == "__main__":
    run()
