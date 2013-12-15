'use strict';

function EditorMainCtrl($scope, $routeParams, $location, Institution) {
    Institution.getModel($routeParams.instID);

    //inherit subscribe from service
    Institution.onDataChange($scope, function(newData) {
        $scope.institution = newData;
    });

    //go to 'edit' View
    $scope.edit = function () {
        $location.path($routeParams.instID + '/edit');
    };

    $scope.setInst = function(instID) {
        // This method is to change without using the location
        // Institution.setInst(instID);
        $location.path(instID);
    }

}


/* Controllers */
function EditorEditCtrl($scope, $routeParams, $location, Institution) {
    Institution.getModel($routeParams.instID);

    Institution.onDataChange($scope, function(obj) {
        $scope.institution = obj;
    });

    //in case we want this to do something more interesting...
    //e.g. copy model here, discard is easier, save is saving the new one
    $scope.$on('formDirty', function (event, formName) {
        console.log("form has been altered");
    });

    $scope.reset = function() {
        Institution.reset();
        //myForm is from name the DOM... may be a better way of handling this (i.e. more dynamic)
        console.log($scope);
        $scope.myForm.$setPristine();
    };

    //save edits, go back to 'view' View
    $scope.save = function() {
        Institution.save($scope.institution);
        $location.path('/' + $routeParams.instID);
    };

    //discard edits, go back to 'view' View
    $scope.discard = function() {
        Institution.discard($scope.institution);
        $location.path('/' + $routeParams.instID);
    };
}