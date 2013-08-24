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
        $scope.editable = {schema: "Schema"};
    };




    //testing - typeahead
    $scope.states = ['Alabama', 'Alaska', 'Arizona', 'Arkansas', 'California', 'Colorado', 'Connecticut', 'Delaware', 'Florida', 'Georgia', 'Hawaii', 'Idaho', 'Illinois', 'Indiana', 'Iowa', 'Kansas', 'Kentucky', 'Louisiana', 'Maine', 'Maryland', 'Massachusetts', 'Michigan', 'Minnesota', 'Mississippi', 'Missouri', 'Montana', 'Nebraska', 'Nevada', 'New Hampshire', 'New Jersey', 'New Mexico', 'New York', 'North Dakota', 'North Carolina', 'Ohio', 'Oklahoma', 'Oregon', 'Pennsylvania', 'Rhode Island', 'South Carolina', 'South Dakota', 'Tennessee', 'Texas', 'Utah', 'Vermont', 'Virginia', 'Washington', 'West Virginia', 'Wisconsin', 'Wyoming'];

}]);
