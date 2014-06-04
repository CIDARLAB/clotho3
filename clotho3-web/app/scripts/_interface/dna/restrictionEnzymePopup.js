angular.module('clotho.clothoDirectives')
/**
 * @ngdoc directive
 * @name restriction-enzyme-popup
 *
 *
 *
 */
	.directive('restrictionEnzymePopup', function ($clothoPopup) {
		return $clothoPopup('restrictionEnzymePopup', {
			placement : 'top',
			trigger : 'mouseenter'
		});
	})
	.directive('restrictionEnzymePopupInner', function () {
		return {
			restrict: 'A',
			replace: true,
			scope: {
				enzyme : '=?sharableModel',
				placement: '@popupPlacement',
				reposition : '&'
			},
			templateUrl: 'views/_interface/dna/restrictionEnzymePopup.html'
		};
	});