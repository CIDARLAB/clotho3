'use strict';

$clotho.extensions.controller('clothoTool_digestCuts', function ($scope, Clotho) {

	Clotho.query({schema: 'org.clothocad.model.RestrictionEnzyme'}, {mute : true}).then(function (results) {
		$scope.enzymes = results;

		//let's activate one to start
		_.find($scope.enzymes, {name : 'EcoRI'})['active'] = true;
	});

	$scope.sequence = 'acaacgtctcacggatccagtcggaattctacatgcatcgatcggctcgactagtgcatgctagctagcggtaccacgcgtgacggatccagatcgactagc';
})
.directive('digestToolHighlight', function (Digest, $parse, $compile) {
	return {
		restrict: 'A',
		scope: {
			enzymes : '=',
			sequence : '='
		},
		link: function (scope, element, attrs, ngModel) {

			function getActiveEnzymes() {
				return _.filter(scope.enzymes, 'active');
			}

			scope.$watch('sequence', function (newval) {
				reprocess();
			});

			scope.$watch(function () {
				return getActiveEnzymes()
			}, function (newval) {
				scope.activeEnzymes = newval;
				reprocess();
			}, true);

			function reprocess() {
				var markedSequence = Digest.markSites(scope.sequence, scope.activeEnzymes);

				//find match sites
				var findMatch = /\((.+?)\)/gi;
				var addedAnnotations = markedSequence.replace(findMatch, function (match, one) {

					var cleanedMatch = one.replace(/[\^_]/g, '');
					var enzymeName = _.find(scope.activeEnzymes, {match : cleanedMatch})['name'];

					return '<digest-tool-match sequence="' + cleanedMatch + '" enzyme-name="' + enzymeName + '">'+one+'</digest-tool-match>'
				});

				console.log(addedAnnotations);

				//replace cut marks
				addedAnnotations = addedAnnotations.replace(/\^/gi, '<digest-tool-cut-top></digest-tool-cut-top>');
				addedAnnotations = addedAnnotations.replace(/_/gi, '<digest-tool-cut-bottom></digest-tool-cut-bottom>');

				element.html($compile('<div>' + addedAnnotations + '</div>')(scope));
			}
		}
	};
})
.directive('digestToolMatch', function () {
	return {
		restrict: 'E',
		transclude: true,
		scope: true,
		template: '<span tooltip="{{ enzymeName }}" tooltip-placement="mouse" tooltip-animation="false" tooltip-append-to-body="true" ng-transclude></span>',
		link : function (scope, element, attrs) {
			element.css({ backgroundColor: '#fcc' });
			scope.enzymeName = attrs.enzymeName;
		}
	};
})
.directive('digestToolCutTop', [function () {
	return {
		restrict: 'E',
		link: function (scope, element, attrs) {
			element.html('&#8595;');
			element.css('color', '#f00');
		}
	};
}])
.directive('digestToolCutBottom', [function () {
	return {
		restrict: 'E',
		link: function (scope, element, attrs) {
			element.html('&#8593;');
			element.css('color', '#00f');
		}
	};
}]);