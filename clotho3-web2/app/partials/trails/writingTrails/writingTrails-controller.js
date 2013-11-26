'use strict';

angular.module('clotho.extensions').controller('writingTrailsCtrl', function($compile, $scope, Clotho, $http, $location) {
    $scope.getSchema = function(schemaName) {
      Clotho.query({schema: "ClothoSchema", name: schemaName})
        .then(function(results) {
          $scope.demoSchema = results[0]
        });
    };

    $scope.getTrail = function(trailName) {
      Clotho.query({schema: 'Trail', name: trailName})
        .then(function(results) {
          $scope.demoTrail = results[0]
        });
    };

    //temporary hack while Trail schema can't be imported
    $scope.getTrailHack = function() {
        $http.get('partials/trails/writingTrails/demo_trailSchema.json').then(function(data) {
            $scope.demoSchema = data.data;
        });
    };

	$scope.$watch(function() {
			return $location.path()
		}, function() {
		$scope.demoSchema = {};
		$scope.demoTrail = {};
	})

})
.controller('writingTrailsUICtrl', function(Clotho, $scope, $focus, $dialog) {
    $scope.videoDemo = {
        id : '<11 char videoID>',
        params : '<youtubeParamsObj or URL>',
        autoplay : 'true|false',
        mini : 'true|false'
    };

    $scope.demoSize = {
        "params" : {
            "height" : 394,
            "width" : 700
        }
    };

    $scope.showSearchTypeOut = function () {
        $focus.typeOutSearch('Demo search text');
    };

    $scope.showDialog = function() {
        $dialog.messageBox('Demo Dialog', 'This is some simple text. It can contain <code>HTML</code>', [{label: "OK", cssClass: "btn-primary", result: true}]).open()
            .then(function() {
                return $dialog.messageBox('Chaining', 'The $dialog and $focus services return promises so you can chain together UI elements based on their completion / response', [{label: "OK", cssClass: "btn-primary", result: true}]).open()
            })
            .then(function() {
                return $focus.typeOutSearch("$focus demo")
            })

    };

    //demo video modal
    $scope.openVideo = function() {
        $dialog.video('ivif214cQLE', {}).open();
    };

    $scope.schemas = [];
    Clotho.query({"schema": "Schema"}).then(function(data) {
        $scope.schemas = data;
    });


    $scope.dropdownItems = [
        "The first choice!",
        "And another choice for you.",
        "but wait! A third!"
    ];


})
.controller('writingTrailsClothoCtrl', function(Clotho, $scope) {

});