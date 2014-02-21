'use strict';

$clotho.extensions.controller('demoClothoDirectives', function($scope, Clotho, $modal) {

    $scope.sequence = 'acgtACGTACTGACacagt';
    $scope.chooseFunction = 'lowercase';

    $scope.clothoFunctions = [];
    Clotho.query({schema: "Function"}).then(function(result) {
        console.log(result);
        $scope.clothoFunctions = result;
    });


    //demo video modal
    $scope.openVideo = function() {
        $modal.video('ivif214cQLE', {});
    }
});