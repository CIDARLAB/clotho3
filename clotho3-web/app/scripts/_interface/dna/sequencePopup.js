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

				function nucseqAsFASTA (nucseq, forceName) {
					return '>' + (forceName || nucseq.name || scope.title) + (nucseq.description ? ' ' + nucseq.description : '') + '\n' + nucseq.sequence;
				}

				scope.$watch('inputModel', function (newval) {
					scope.model = angular.isArray(newval) ? newval : [newval];
					if (newval) {
						var FASTAs = angular.map(scope.model, function (nucseq, index) {
							var name = nucseq.name || ( (scope.title ? scope.title + ' - ' : '') + 'Fragment ' + (index + 1) );
							return nucseqAsFASTA(nucseq, name);
						});
						scope.generatedFASTA = FASTAs.join('\n\n');
					}
				});
			}
		};
	});