import graphene
import viewer.schema, hypothesis.schema, hotspots.schema, pandda.schema, scoring.schema


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
