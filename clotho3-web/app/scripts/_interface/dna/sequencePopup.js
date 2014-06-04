angular.module('clotho.clothoDirectives')
/**
 * @ngdoc directive
 * @name sequence-popup
 *
 * @description can handle a single NucSeq or array of NucSeqs
 *
 */
	.directive('sequencePopup', function ($clothoPopup) {
		return $clothoPopup('sequencePopup', {
			placement : 'top',
			trigger : 'mouseenter'
		});
	})
	.directive('sequencePopupInner', function () {
		return {
			restrict: 'A',
			replace: true,
			scope: {
				inputModel : '=?sharableModel',
				title : '=popupTitle',
				placement: '@popupPlacement',
				reposition : '&'
			},
			templateUrl: 'views/_interface/dna/sequencePopup.html',
			link: function (scope, element, attrs) {
				scope.$watch('inputModel', function (newval) {
					scope.model = angular.isArray(newval) ? newval : [newval];
				});
			}
		};
	});