import click
import os
import json
import pathlib
from typing import Optional
from tqdm import tqdm
import openfe

from loguru import logger
from gufe import tokenization
from alchemiscale import AlchemiscaleClient, Scope, ScopedKey


def get_network(filename: pathlib.Path) -> openfe.AlchemicalNetwork:
    """
    Read a serialized AlchemicalNetwork

    Parameters
    ----------
    filename : pathlib.Path
      A path to a file with the serialized AlchemicalNetwork

    Returns
    -------
    AlchemicalNetwork
      An AlchemicalNetwork with the desired Transformations.
    """
    return openfe.AlchemicalNetwork.from_json(filename)


@click.command
@click.option(
    '--network_filename',
    type=click.Path(dir_okay=False, file_okay=True, path_type=pathlib.Path),
    required=True,
    help="Path to the serialized AlchemicalNetwork",
)
@click.option(
    '--org_scope',
    type=str,
    required=True,
    help='The organization scope name',
)
@click.option(
    '--scope_name_campaign',
    type=str,
    required=True,
    help='The campaign transformation scope name',
)
@click.option(
    '--scope_name_project',
    type=str,
    required=True,
    help='The project transformation scope name',
)
@click.option(
   '--repeats',
   type=int,
   required=True,
   help='The number of repeats per transformation',
)
@click.option(
    '--scopekey_output',
    type=click.Path(dir_okay=False, file_okay=True, path_type=pathlib.Path),
    required=False,
    default="scoped-key.dat",
    help="The file name for where we write the scope key",
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
    default=0.0001,
    help="Set the priority of the network (network weight)",
)
def run(
    network_filename: pathlib.Path,
    org_scope: str,
    scope_name_campaign: str,
    scope_name_project: str,
    repeats: int,
    scopekey_output: pathlib.Path,
    user_id: Optional[str],
    user_key: Optional[str],
    priority: float,
):
    """
    Submit and action an AlchemicalNetwork on Alchemiscale.

    Parameters
    ----------
    network_filename : pathlib.Path
      Path to the serialized AlchemicalNetwork.
    org_scope : str
      The organization Scope name.
    scope_name_campaign : str
      The campaign Scope name.
    scope_name_project : str
      The project Scope name.
    repeats : int
      The number of repeats to action per task.
    scopekey_output : pathlib.Path
      A path to where to write the serialized ScopeKey.
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

    # Get the alchemical network
    alchemical_network = get_network(network_filename)
    logger.info(alchemical_network)

    # Set the scope for the transformation
    scope = Scope(org_scope, scope_name_campaign, scope_name_project)
    logger.info(f"Scope is set to: {scope}")

    # Create a network and get a scope key
    an_sk = asc.create_network(alchemical_network, scope)

    # When a compute node is ready for a task, it first stochastically
    # chooses a network, then stochastically chooses a task within that network
    # The weight is the likelihood that that network/task will be chosen
    # Choosing a small value allows other networks to take priority on shared
    # compute nodes.
    asc.set_network_weight(an_sk, priority)

    # store the scoped key
    with open(scopekey_output, 'w') as f:
        f.write(str(an_sk))

    # action out tasks
    logger.info(f"Actioning {repeats} repeats for {len(alchemical_network.edges)} transformation")
    for transform_sk in tqdm(asc.get_network_transformations(an_sk)):
        tasks = asc.create_tasks(transform_sk, count=repeats) # will always create a new task
        asc.action_tasks(tasks, an_sk) # idempotent

    # check what the network status looks like
    asc.get_network_status(an_sk)


if __name__ == "__main__":
    run()
