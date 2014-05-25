'use strict';

$clotho.extensions.controller('constructionFiles_wetLabSyntaxCtrl', function($scope, $http, Clotho) {

	//todo - add raw versions

	Clotho.query({schema : "org.clothocad.model.ConstructionFileRaw"}, {mute : true}).then(function(data) {
		$scope.constructionFiles = data
	});

});