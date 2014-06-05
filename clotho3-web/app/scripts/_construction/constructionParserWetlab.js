angular.module('clotho.construction')
/**
 * @ngdoc directive
 * @name constructionParserWetlab
 *
 * @restrict A
 *
 * @description
 * Template a construction file to a wet lab protocol
 *
 * @example
 *
 */
	.directive('constructionParserWetlab', function () {
		return {
			restrict : 'A',
			scope: {
				file : '=constructionParserWetlab'
			},
			templateUrl : 'views/_construction/constructionParserWetlab.html',
			link : function (scope, element, attrs) {}
		}
	})
	.directive('constructionParserWetlabStep', function (ConstructionReactions, $compile) {
		return {
			restrict : 'A',
			link: function (scope, element, attrs) {
				attrs.$observe('stepType', function (newval) {
					element.html($compile('' +
							'<td>' + ConstructionReactions[newval].template_wetlabL + '</td>' +
							'<td>' + ConstructionReactions[newval].template_wetlabR + '</td>'
					)(scope));
				})
			}
		}
	});