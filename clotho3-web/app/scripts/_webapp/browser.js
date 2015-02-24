'use strict';

angular.module('clotho.webapp').controller('BrowserCtrl', function($scope, Clotho, $filter, ClothoSchemas) {

    //todo - pending recent() (see GH#410), incorporate recent list of items

    /* ordering */

    $scope.orderers = [
        {
            name : "Name",
            criteria : "name",
            class: "glyphicon-sort-by-alphabet"
        },
        {
            name : "type",
            criteria : ClothoSchemas.dirtyDetermineType,
            class: "glyphicon-th-list"
        },
        {
            name : "Schema",
            criteria : "schema",
            class : "glyphicon-cog"
        }
    ];

    $scope.setOrder = function (orderer) {
        $scope.currentOrder = (orderer === $scope.currentOrder) ? null : orderer;
    };

    /* filters */

    $scope.filters = [
        {
            name : "Made by me",
            filter : function (sharable) {
                //todo - update pending users
                return angular.isDefined(sharable.author) && sharable.author == true;
            },
            "class" : "glyphicon glyphicon-user"
        },
        {
            name : "Has Description",
            filter : function (sharable) {
                return angular.isDefined(sharable.description) && sharable.description;
            },
            "class" : "glyphicon glyphicon-comment"
        }
    ];

    $scope.currentFilter = function () {return true};
    $scope.setFilter = function (filter) {
        $scope.currentFilter = (filter === $scope.currentFilter) ?
            function () {return true} :
            filter;
    };



    //todo - update when pulling from server
    $scope.collections = [
        {
            name : "My Collection",
            author : "uniqueUserID",
            description : "my Collection of all the things.",
            items : {
                "clotho.developer.maxbates" : {
                    "note1" : "blah blah blah blah"
                },
                "clotho.enzyme.BglII" : {
                    "note1" : "yad yad ayayayayy"
                },
                "clotho.part.jtk2134" : {
                    "note2" : "bling bling blang"
                }
            }
        }
    ];

    $scope.setCollection = function (coll) {
        if ($scope.currentQuery == coll) {
            return;
        }
        $scope.currentQuery = coll;
        $scope.resultArray = [];

        angular.forEach(coll.items, function (notes, id) {
         Clotho.get(id, {mute : true}).then(function (result) {
             $scope.resultArray.push(result);
         });
        });
    };

    //todo - store on server and pull
    $scope.queries = [
        {
            name: "Random",
            query : {}
        },
        {
            name : "Schemas",
            query : {
                schema: "org.clothocad.core.schema.Schema"
            }
        },
        {
            name: "NucSeqs",
            query : {
                "schema" : "org.clothocad.model.NucSeq"
            }
        },
        {
            name: "Vectors",
            query : {
                "schema" : "Vector"
            }
        },
        {
            name : "Contains pBAC",
            query : {
                name : {"$regex"  : 'pBAC'}
            }
        }
    ];

    $scope.newQuery = {};
    $scope.saveNewQuery = function () {
        //todo - persist on server (requires users)
        $scope.queries.push($scope.newQuery);
        $scope.newQuery = {};
    };

    /* query construction */

    $scope.setCurrentQuery = function(value, limit) {
        if ($scope.currentQuery == value) {
            return;
        }
        limit = limit || 200;
        $scope.currentQuery = value;
        Clotho.query(value, {maxResults : limit}).then(function (result) {
            $scope.resultArray = result;
        });
    };

    $scope.sort = function(byCat) {
        if (byCat) {
            if (!!$scope.catSort) {return;}

            $scope.catSort = true;
            $scope.recent = $filter('categorize')($scope.recent_array, 'type');
            $scope.recent['Instance'] = $filter('categorize')($scope.recent['Instance'], 'schema.name');
        } else {
            if (!$scope.catSort && angular.isDefined($scope.catSort)) {return;}

            $scope.catSort = false;
            $scope.recent = {
                "entries" : angular.copy($scope.recent_array)
            }
        }
    };

    /* styling */

    $scope.collectionIconClass = 'glyphicon glyphicon-briefcase';
    $scope.filterIconClass = 'glyphicon glyphicon-filter';
    $scope.queryIconClass = 'glyphicon glyphicon-wrench';

    //init

    $scope.setCurrentQuery($scope.queries[0].query);

});
