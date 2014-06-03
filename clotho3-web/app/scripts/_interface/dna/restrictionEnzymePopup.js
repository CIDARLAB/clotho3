angular.module('clotho.clothoDirectives')
/**
 * @ngdoc directive
 * @name restriction-enzyme-popup
 *
 *
 *
 */
	.directive('restrictionEnzymePopup', function ($clothoPopup) {
		return $clothoPopup('restrictionEnzymePopup');
	})
	.directive('restrictionEnzymePopupInner', function () {
		return {
			restrict: 'A',
			replace: true,
			scope: {
				enzyme : '=?sharableModel',
				placement: '@',
				reposition : '&'
			},
			templateUrl: 'views/_interface/dna/restrictionEnzymePopup.html'
		};
	});