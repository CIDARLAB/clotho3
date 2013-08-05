'use strict';

Application.Editor.controller('EditorCtrl', ['$scope', '$routeParams', '$location', 'Clotho', function($scope, $routeParams, $location, Clotho) {

    //fixme - change path, don't reinstantiate controllers etc.
    //note - could use search param, and set reloadOnSearch to false
    /*$scope.$watch('id', function (newval, oldval) {
       $location.url('/editor/' + $scope.id).replace();
    });*/

    //init()
    $scope.id = $routeParams.id;

    $scope.editableList = [
        {"name" : "Person", "type" : "Sharable", "value" : "sharable_person"},
        {"name" : "Institution", "type" : "Sharable", "value" : "sharable_institution"},
        {"name" : "First Function", "type" : "Function", "value" : "func_first"},
        {"name" : "Second Function", "type" : "Function", "value" : "func_second"}
    ];


//testing - typeahead
    $scope.selected = undefined;
    $scope.states = ['Alabama', 'Alaska', 'Arizona', 'Arkansas', 'California', 'Colorado', 'Connecticut', 'Delaware', 'Florida', 'Georgia', 'Hawaii', 'Idaho', 'Illinois', 'Indiana', 'Iowa', 'Kansas', 'Kentucky', 'Louisiana', 'Maine', 'Maryland', 'Massachusetts', 'Michigan', 'Minnesota', 'Mississippi', 'Missouri', 'Montana', 'Nebraska', 'Nevada', 'New Hampshire', 'New Jersey', 'New Mexico', 'New York', 'North Dakota', 'North Carolina', 'Ohio', 'Oklahoma', 'Oregon', 'Pennsylvania', 'Rhode Island', 'South Carolina', 'South Dakota', 'Tennessee', 'Texas', 'Utah', 'Vermont', 'Virginia', 'Washington', 'West Virginia', 'Wisconsin', 'Wyoming'];



}]);