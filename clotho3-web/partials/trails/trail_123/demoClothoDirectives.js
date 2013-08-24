'use strict';

Application.Extensions.controller('demoClothoDirectives', ['$scope', 'Clotho', '$dialog', function($scope, Clotho, $dialog) {

    $scope.sequence = 'acgtACGTACTGACacagt';
    $scope.chooseFunction = 'lowercase';

    $scope.clothoFunctions = [];
    Clotho.query({schema: "Function"}).then(function(result) {
        console.log(result);
        $scope.clothoFunctions = result;
    });


    //demo video modal
    $scope.openVideo = function() {
        $dialog.video('ivif214cQLE', {}).open();
    }
}]);