'use strict';

$clotho.extensions.controller('clothoIntro_SimulationCtrl', ['$scope', '$focus', '$timeout', '$dialog', '$http', 'DNA', 'Digest', 'PCR', 'Construction', 'Clotho', function($scope, $focus, $timeout, $dialog, $http, DNA, Digest, PCR, Construction, Clotho) {

    $http.get('models/construction_gfp.json').then(function(data) {
        $scope.constructionFile = data.data
    });

    $scope.reactions = Construction.reactions;

}]);