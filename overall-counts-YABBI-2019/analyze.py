"""
dsw7@sfu.ca
The analysis script that processes the results of the Met-aromatic algorithm.
I used this script to process:

    -- ma.non_redundant_no_ang_limit
    -- ma.non_redundant_1095_ang_limit

And generate:

    -- ma.bridges_1095_ang_limit
    -- ma.bridges_no_ang_limit
    -- ma.non_redundant_1095_ang_limit
    -- ma.non_redundant_no_ang_limit
    -- ma.pairs_1095_ang_limit
    -- ma.pairs_no_ang_limit

"""

# -------------------------------------------------------------------------------------------
from pymongo import MongoClient, errors
from networkx import Graph, connected_components
from itertools import groupby

# -------------------------------------------------------------------------------------------
# manually input database and collection names

database = "ma"
collection = "non_redundant_1095_ang_limit"  # "non_redundant_no_ang_limit"
collection_pairs = "pairs_1095_ang_limit"  # "pairs_no_ang_limit"
collection_bridges = "bridges_1095_ang_limit"  # "bridges_no_ang_limit"


# -------------------------------------------------------------------------------------------
mongoport = 27017
mongohost = "localhost"
client = MongoClient(mongohost, mongoport)


# -------------------------------------------------------------------------------------------


def count_entries(query_database, query_collection):
    """
    :param query_database: Database in which Met-aromatic pairs are located.
    :param query_collection: Collection in which Met-aromatic pairs are located.
    :return: Nothing. Function performs operation pass-by style.
    """
    cursor = client[query_database][query_collection]
    count = len(cursor.distinct('code'))
    if not count:
        exit(' -- Empty collection: {}. Exiting.'.format(query_collection))
    else:
        print(' -- Number of unique entries analyzed: {}\n'.format(count))


def breakdowns_by_order(query_database, query_collection):
    """
    :param query_database: Database in which Met-aromatic pairs are located.
    :param query_collection: Collection in which Met-aromatic pairs are located.
    :return: Nothing. Function performs operation pass-by style.
    """
    query = [
        {
            '$group': {
                '_id': {
                    'code': "$code",
                    'arores': "$arores",
                    'met': "$met"
                },
                'order': {
                    '$sum': 1
                }
            }
        },
        {
            '$project': {
                'order': 1
            }
        },
        {
            '$group': {
                '_id': '$order',
                'count': {
                    '$sum': 1
                },
            }
        },
        {
            '$sort': {
                '_id': 1
            }
        }
    ]

    print(' -- Breakdowns by order: ')
    cursor = client[query_database][query_collection].aggregate(query)
    for entry in cursor:
        print(' -- Order: {} | Count: {}'.format(entry.get('_id'), entry.get('count')))


def get_pairs(query_database, query_collection, name_collection_outgoing):
    """
    Function gets methionine-aromatic pairs from base collection obtained from
    runner.py. Function then exports the pairs to a new collection.

    :param query_database: Database in which Met-aromatic pairs are located.
    :param query_collection: Collection in which Met-aromatic pairs are located.
    :param name_collection_outgoing: The name of the MongoDB collection to export to.
    :return: Nothing. Function performs operation pass-by style.
    """
    # ensure duplicates are not being loaded
    client[query_database][name_collection_outgoing].drop()

    query = [
        {
            '$group': {
                '_id': {
                    'code': "$code",
                    'EC': {'$substr': ["$ec", 0, 1]}  # 1.2.3.56 -> 1
                },
                'pairs': {
                    '$addToSet':  {
                        '$concat': ["$aro", "", "$arores", "|", "MET", "", "$met"]
                    }
                }
            }
        },
        {'$project': {'pairs': 1, 'EC': '$_id.EC', 'code': '$_id.code', '_id': 0}},
        {'$out': name_collection_outgoing}
    ]

    client[query_database][query_collection].aggregate(query)
    print("\n -- Wrote Met-aromatic pair data to collection: {}".format(name_collection_outgoing))


def get_bridges_from_pairs(query_database, query_collection, name_collection_outgoing, n=2):
    """
    Function gets pairs from collection generated in get_pairs function. Function
    then takes advantage of NetworkX to find n-bridges using network theory approaches.
    Function then exports to another collection in MongoDB: a bridges collection.

    :param query_database: Database in which Met-aromatic pairs are located.
    :param query_collection: Collection in which Met-aromatic pairs are located.
    :param name_collection_outgoing: The MongoDB collection to export data to.
    :param n: 2-bridge, 3-bridge, 4-bridge, ..., n-bridge
    :return: Nothing. Function performs operation pass-by style.
    """
    if n < 2:
        exit("Incorrect bridge order. A bridge must be of n >= 2!")
    else:
        # ensure duplicates are not being loaded
        client[query_database][name_collection_outgoing].drop()

        bridges = []
        for entry in client[query_database][query_collection].find():
            pairs = []
            for pair in entry.get('pairs'):
                pairs.append(tuple(pair.split('|')))  # 'TYR123|MET123' -> ('TYR123', 'MET123')

            G = Graph()
            G.add_edges_from(pairs)

            for disconnects in list(connected_components(G)):
                if len(disconnects) == n + 1:
                    if ''.join(disconnects).count('MET') == 1:  # remove inverse bridges -> MET-ARO-MET
                        bridges.append(
                            {
                                'code': entry.get('code'),
                                'EC': entry.get('EC'),
                                'bridge': list(disconnects)
                            }
                        )
                    else:
                        pass
                else:
                    pass

        try:
            client[query_database][name_collection_outgoing].insert_many(bridges)
        except errors.BulkWriteError as pymongo_exception:
            print(pymongo_exception.details['writeErrors'])
        else:
            print(' -- Collected all {}-bridges and exported to collection: {} \n'.format(n, name_collection_outgoing))


def count_bridges(query_database, query_collection):
    """
    Function gets bridges from MongoDB bridges collection generated by
    get_bridges_from_pairs function and performs a count.

    :param query_database: Database in which bridges are located.
    :param query_collection: Collection in which bridges are located.
    :return: Count of bridges by type. A dict of form:
    {
        PHE-PHE: 3,
        TYR-TYR: 0,
        TRP-TRP: 0,
        TYR-TRP: 0,
        PHE-TYR: 2,
        PHE-TRP: 1
    }
    """
    bridges = [i.get('bridge') for i in client[query_database][query_collection].find()]

    # clean up data prior to count
    refined = []
    for bridge in bridges:
        bridge.remove(*(amino_acid for amino_acid in bridge if "MET" in amino_acid))
        refined.append(''.join([s for s in '-'.join(bridge) if not s.isdigit()]))

    # actual counts
    FF = refined.count('PHE-PHE')
    YY = refined.count('TYR-TYR')
    WW = refined.count('TRP-TRP')
    YW = refined.count('TYR-TRP') + refined.count('TRP-TYR')
    FY = refined.count('TYR-PHE') + refined.count('PHE-TYR')
    FW = refined.count('PHE-TRP') + refined.count('TRP-PHE')

    return {
        'PHE-PHE': FF,
        'TYR-TYR': YY,
        'TRP-TRP': WW,
        'TYR-TRP': YW,
        'PHE-TYR': FY,
        'PHE-TRP': FW
    }


def count_bridges_by_ec(query_database, query_collection):
    """
    Group by EC classifier and then count bridges. Custom Jeff Warren request.

    :param query_database: Database in which bridges are located.
    :param query_collection: Collection in which bridges are located.
    :return: A dict of dicts for containing bridge counts grouped by EC codes.
    """
    bridges = [row for row in client[query_database][query_collection].find()]

    # need to sort prior to grouping with itertools.groupby
    key = lambda column: column['EC']
    bridges.sort(key=key)

    results = {}

    for ec, bridge in groupby(bridges, key=key):
        ec = '0' if ec == '' else ec  # replace '' with '0'
        row, refined = [i.get('bridge') for i in list(bridge)], []
        for r in row:
            r.remove(*(amino_acid for amino_acid in r if "MET" in amino_acid))
            refined.append(''.join([s for s in '-'.join(r) if not s.isdigit()]))

        # actual counts
        FF = refined.count('PHE-PHE')
        YY = refined.count('TYR-TYR')
        WW = refined.count('TRP-TRP')
        YW = refined.count('TYR-TRP') + refined.count('TRP-TYR')
        FY = refined.count('TYR-PHE') + refined.count('PHE-TYR')
        FW = refined.count('PHE-TRP') + refined.count('TRP-PHE')

        results[ec] = {
            'PHE-PHE': FF,
            'TYR-TYR': YY,
            'TRP-TRP': WW,
            'TYR-TRP': YW,
            'PHE-TYR': FY,
            'PHE-TRP': FW
        }

    return results


# -------------------------------------------------------------------------------------------


if __name__ == '__main__':
    print(' =============')
    print(' --- START ---')
    print(' =============')
    print('\n')

    print(' -- Connected to MongoDB on:')
    print(' -- Port: {}'.format(mongoport))
    print(' -- Host: {}'.format(mongohost))
    print(' -- Database: {}'.format(database))
    print(' -- Base collection: {}'.format(collection))
    print(' -- Met-aromatic pairs collection: {}'.format(collection_pairs))
    print(' -- Bridges collection: {}\n'.format(collection_bridges))
    print(' -- Analyzing...\n')

    count_entries(database, collection)
    breakdowns_by_order(database, collection)
    get_pairs(database, collection, collection_pairs)
    get_bridges_from_pairs(database, collection_pairs, collection_bridges)

    print(" -- Bridges: ")
    bridges_overall = count_bridges(database, collection_bridges)
    for k in bridges_overall:
        print(' -- {} | {}'.format(k, bridges_overall.get(k)))

    print("\n -- Bridge counts by EC: ")
    print(" -- EC | { counts }")
    bridges_by_EC = count_bridges_by_ec(database, collection_bridges)
    for k in bridges_by_EC:
        print(' --  {} | {}'.format(k, bridges_by_EC.get(k)))

    print('\n')
    print(' =============')
    print(' ---- END ----')
    print(' =============')

client.close()
