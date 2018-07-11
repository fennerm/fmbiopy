import click
from plumbum import local


class PlumbumPath(click.Path):
    """A click parameter type which auto converts to plumbum LocalPath."""

    name = "plumbum_path"

    def convert(self, value, param, ctx):
        return local.path(super().convert(value, param, ctx))
