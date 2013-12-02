'use strict';

Application.Editor.controller('EditorCtrl', ['$scope', '$routeParams', '$location', 'Clotho', function($scope, $routeParams, $location, Clotho) {

    //fixme - change path, don't reinstantiate controllers etc.
    //note - could use search param, and set reloadOnSearch to false
    /*$scope.$watch('id', function (newval, oldval) {
       $location.url('/editor/' + $scope.id).replace();
    });*/

    //init()
    $scope.id = $routeParams.id;

    //data

    $scope.schemas = [];
    Clotho.query({"schema": "Schema"}).then(function(data) {
        $scope.schemas = data;
    });

    $scope.editableList = [];
    Clotho.query({}).then(function (data){
        var editableList = [];
        for (var i=0; i<data.length; i++){
            if (data[i].schema != "BuiltInSchema"){
                editableList.push(data[i]);
            }
        }
        $scope.editableList = editableList;
    });

    //functions

    $scope.createNewObject = function() {
        Clotho.create({schema: $scope.selected}).then(function(id) {
            console.log(id);
            $scope.id = id;
        });
        $scope.selected = undefined;
    };

    $scope.logSelected = function() {
        console.log($scope.selected);
    };

    $scope.createNewSchema = function() {
        $scope.editable = {schema: "ClothoSchema", language: "JSONSCHEMA"};
    };

}]);