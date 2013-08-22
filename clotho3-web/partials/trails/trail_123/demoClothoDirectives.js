'use strict';

Application.Extensions.controller('demoClothoDirectives', ['$scope', 'Clotho', function($scope, Clotho) {

    $scope.sequence = 'acgtACGTACTGACacagt';
    $scope.chooseFunction = 'lowercase';

    $scope.clothoFunctions = [];
    Clotho.query({schema: "Function"}).then(function(result) {
        console.log(result);
        $scope.clothoFunctions = result;
    });
}]);