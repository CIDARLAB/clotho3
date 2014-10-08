'use strict';

/**
 *
 * @attr clothoSharable {Sharable} Sharable to show info for
 */

angular.module('clotho.interface')
	.directive('clothoSharable', function ($injector, Clotho, ClothoSchemas) {

		function wrapInnerUrl(type) {
			return 'views/_interface/sharables/' + (type || 'default') + '.html';
		}

		var editorPresent = $injector.has('clothoEditorDirective');

		return {
			restrict: "A",
			scope: {
				sharable: '=clothoSharable'
			},
			templateUrl: 'views/_interface/clothoSharable.html',
			link: function clothoSharableLink(scope, element, attrs, controller) {

				//future - DRY: a lot of this is duplicated between this + sharablePopup
				scope.$watch('sharable', function (newval) {
					//note - dirty check type here to minimize requests, pending GH#418 and GH#416
					scope.type = ClothoSchemas.dirtyDetermineType(newval);
					scope.iconClass = ClothoSchemas.determineSharableIcon(scope.type);
					scope.labelClass = 'label-' + ClothoSchemas.typeToColorClass(scope.type);
					scope.schemaName = ClothoSchemas.mapSchemaIdToName(newval.schema);
				});


				scope.editorPresent = editorPresent;
				scope.edit = Clotho.edit;


				/*
				//todo - use ClothoSchemas to determine type, use specific template
				//stale!
				ClothoSchemas.determineSharableType(scope.sharable.id).then(function (detType) {
					$http.get(wrapInnerUrl(detType), {cache: $templateCache})
						.success(function (data, status, headers, config) {
							element.html(data);
							$compile(element.contents())(scope);
						})
						.error(function (data, status, headers, config) {
							$http.get(wrapInnerUrl(), {cache: $templateCache})
								.success(function (data, status, headers, config) {
									//force storage for template that didn't exist
									$templateCache.put(wrapInnerUrl, data);

									element.html(data);
									$compile(element.contents())(scope);
								});
						});
					});
					*/
			}
		};
	});