import graphene

import hotspots.schema
import hypothesis.schema
import pandda.schema
import scoring.schema
import viewer.schema


class Query(
    viewer.schema.Query,
    hypothesis.schema.Query,
    hotspots.schema.Query,
    pandda.schema.Query,
    scoring.schema.Query,
    graphene.ObjectType,
):
    # This class will inherit from multiple Queries
    # as we begin to add more apps to our project
    pass


schema = graphene.Schema(query=Query)
