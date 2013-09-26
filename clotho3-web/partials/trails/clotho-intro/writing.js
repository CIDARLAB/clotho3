'use strict';

Application.Extensions.controller('clothoIntro_WritingCtrl', ['$scope', '$focus', '$timeout', '$dialog', '$http', 'DNA', 'Digest', 'PCR', 'Construction', 'Clotho', function($scope, $focus, $timeout, $dialog, $http, DNA, Digest, PCR, Construction, Clotho) {

    $scope.getConstructionFile = function(selection) {

        var files = [
            '/models/construction_demo.json',
            '/models/construction_gfp.json',
            '/models/construction_vio.json'
        ];

        $http.get(files[selection]).then(function(data) {
            $scope.constructionFile = data.data
        });
    };

    $scope.getConstructionFile(2);

    $scope.reactions = Construction.reactions;

}]);