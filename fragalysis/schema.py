import graphene
import pandda.schema

import hotspots.schema
import hypothesis.schema
import scoring.schema
import viewer.schema
import xcdb.schema


class Query(
    viewer.schema.Query,
    hypothesis.schema.Query,
    hotspots.schema.Query,
    pandda.schema.Query,
    scoring.schema.Query,
    xcdb.schema.Query,
    graphene.ObjectType,
):
    # This class will inherit from multiple Queries
    # as we begin to add more apps to our project
    pass


schema = graphene.Schema(query=Query)
