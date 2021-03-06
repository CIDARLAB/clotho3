'use strict';

angular.module('clotho.webapp').controller('EditorCtrl', function ($scope, $route, $location, Clotho, ClothoSchemas) {

    /* route handling */

    //check for updates to editable's id, update URL
    $scope.$watch('editable.id', function (newval, oldval) {
        $location.search('id', newval || null).replace();
    });

    //check for updates to URL - this will only be triggered when URL is not set by above action
    $scope.$on('$routeUpdate', function(scope, next, current) {
        var updateId = next.params.id;

        // only activate if new -- won't happen if dropdown changed this already
        if (!angular.isEmpty($scope.editable) && $scope.editable.id != updateId) {
            $scope.editableId = updateId;
        }
    });

    $scope.$on('$destroy', function () {
        $location.search('id', null).replace();
    });

    //initial query handling
    if ( $route.current.params.id ) {
        $scope.editableId = $route.current.params.id;
    }

    /*
    initial route handling:

        /editor/query/:queryTerm

    currently not supporting
    todo - remove this part of route when change search term
    */
    /*
    var queryResult = $route.current.locals.queryResult;
    queryResult && console.log('query result', queryResult);
    if (angular.isDefined(queryResult) && queryResult.length) {
        $scope.editable = queryResult[0];
        console.log('\n\n\nEDITOR CONTROLLER', $scope.editable);
    }
  */

    $scope.editModePass = false;

    //data

    $scope.objectTypes = ClothoSchemas.sharableTypes;

    $scope.schemas = [];
    ClothoSchemas.retrievedSchemas.then(function (schemas) {
        $scope.schemas = schemas;
    });

    //functionality

    $scope.queryWrapper = function (text) {
        return Clotho.autocomplete(text).then(function (results) {
            return results || [];
        });
    };

    $scope.createNewNonInstance = function (type) {
        $scope.editable = ClothoSchemas.createScaffold(type);
        $scope.editModePass = true;
    };

    $scope.createNewInstance = function (item, model, label) {
        $scope.editable = {"id": item.id};
        $scope.editModePass = true;
    };

    $scope.editExisting = function (item, query) {
        $scope.editableId = item.id;
        $scope.editModePass = true;
    };
});
