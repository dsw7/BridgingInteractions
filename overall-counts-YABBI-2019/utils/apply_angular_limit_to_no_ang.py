"""
dsw7@sfu.ca
A quick script I used for applying an angular limit to the "non_redundant_no_ang_limit" collection.
This was a JJW request. Generates the new collection:

    -- ma.non_redundant_1095_ang_limit

"""

# -------------------------------------------------------------------------------------------
from pymongo import MongoClient

mongoport = 27017
mongohost = "localhost"
database = "ma"
client = MongoClient(mongohost, mongoport)


def apply_angular_condition(query_database, name_collection_out="non_redundant_1095_ang_limit", cutoff=109.5):
    """
    Function for applying angular condition to non_redundant_no_ang_limit collection. Recall that no angular
    condition was applied to this collection during data mining.

    :param query_database: Database in which Met-aromatic pairs are located
    :param name_collection_out: The MongoDB collection to export data to
    :param cutoff: Angular cutoff to apply to Met-theta / Met-phi
    :return: A new collection where "result" field is a boolean indicating whether condition is true or false
    """

    query = [
        {
            "$match": {
                "$or": [
                    {
                        "met-theta": {"$lte": cutoff}
                    },
                    {
                        "met-phi": {"$lte": cutoff}
                    },
                ]
            }
        },
        {
            "$out": name_collection_out
        }
    ]

    client[query_database]['non_redundant_no_ang_limit'].aggregate(query)


if __name__ == '__main__':
    apply_angular_condition(database)

client.close()
