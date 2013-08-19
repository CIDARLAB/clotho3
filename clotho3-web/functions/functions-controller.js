'use strict';

Application.Functions.controller('FunctionsCtrl', ['$scope', 'Clotho', '$filter', function($scope, Clotho, $filter) {

//init

    $scope.nav = ["Browse", "Create"];

    $scope.setCurrent = function(mode) {
        $scope.current = mode;
    };


    $scope.newDescription = function(){
        return  {
            name: "NewFunction",
            description: "",
            language: "JavaScript",
            source: "",
            canDoIt: ""
        };
    };

    $scope.submit = function (description){
        //reference
        
        description.schema = "JavaScriptFunction";
        
        Clotho.create(description).then(function (data){
            Clotho.get(data).then(function (data){
                $scope.schemas.push(data);
                $scope.created = $scope.newDescription();
            })
        });


    };


    $scope.LANGUAGES = ["JavaScript", "Python", "Java"];
        Clotho.query({"schema":"Function"}).then(function(data){
            $scope.functions = data;
        });
    
    $scope.setCurrent($scope.nav[0]);
    $scope.created = $scope.newDescription();
}]);

