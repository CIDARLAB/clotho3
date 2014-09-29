angular.module('clotho.trails').service('Trails', function(Clotho, $q, $location) {

  // download dependencies object - css, mixin, script - return function to for onload
  var downloadDependencies = $clotho.extensions.downloadDependencies;

  var compile = function TrailCompile(trail) {

		//If pass by reference (depending on Clotho.get() ) need to copy to don't edit in dependencies
		//trail = angular.copy(trail);

		var transcludes = trail.dependencies || null,
			deferred = $q.defer();

		if (!transcludes && !trail.mixin) {
			deferred.resolve(trail);
			return deferred.promise;
		}

		var final_contents = [],
			promises = [];

		//get the transcluded trails... will be fast if in collector already
		(transcludes) && angular.forEach(transcludes, function(id) {
			promises.push(Clotho.get(id));
		});

		//download dependencies and transcluded trails
		downloadDependencies(trail.dependencies)
			.then(function () {
				//keep separate so the next step in the chain has the trails to transclude
				return $q.all(promises);
			})
			//after download all, pluck out the chapters we need
			.then(function (downloads) {

				//reorganize transcludes so can reference by id
				transcludes = {};
				angular.forEach(downloads, function(transclude) {
					transcludes[transclude.id] = transclude;
				});

				//iterate through trail, pushing in chapters
				angular.forEach(trail.contents, function (chapter, ind) {

					if (angular.isUndefined(chapter.transclude)) {
						final_contents.push(chapter);
					} else {
						//chapters to include :
						var chapterID = chapter.transclude.id,
							chapterNum = chapter.transclude.chapters;

						if ((chapterNum == "all") || (typeof chapterNum == 'undefined')) {
							for (var i = 0; i < transcludes[chapterID]['contents'].length; i++) {
								final_contents.push(transcludes[chapterID]['contents'][i]);
							}
						} else {
							var startStop = chapterNum.split("-");
							if (startStop.length == 1) {
								final_contents.push(transcludes[chapterID]['contents'][startStop[0]]);
							} else {
								if (startStop[0] > startStop[1])
									return "wrong format - start must be smaller than end";

								for (var i = startStop[0]; i <= startStop[1]; i++) {
									final_contents.push(transcludes[chapterID]['contents'][i])
								}
							}
						}
					}
				});

				trail.contents = final_contents;
				deferred.resolve(trail);
			});

		return deferred.promise;
	};

	//in form <Chapter>-<page>
	var pageExists = function trailPageExists (Trail, indices) {
		if (angular.isUndefined(indices)) {
			return false;
		}
		var pos = angular.isString(indices) ? indices.split("-") : [0,0];
		var chapter = Trail.contents[pos[0]];
		if (angular.isUndefined(chapter)) {
			return false;
		}
		var page = chapter['pages'][pos[1]];
		return angular.isDefined(page);
	};

	//in form <Chapter>-<Page>
	var extractPage = function trailExtractPage (Trail, indices) {
		var pos = angular.isString(indices) ? indices.split("-") : [0,0];
		if (pageExists(Trail, indices)) {
			return Trail.contents[pos[0]]['pages'][pos[1]];
		} else {
			return null;
		}
	};

	//bring back in logic from trail-module.js
	var calcNextPage = function trailNextPage(Trail, oldpos) {
		oldpos = angular.isString(oldpos) ? oldpos.split("-") : [0, -1];
		var newpos;

		//check next page
		if (typeof Trail.contents[oldpos[0]]['pages'][+oldpos[1] + 1] != 'undefined')
			newpos = oldpos[0] + '-' + (+oldpos[1] + 1);
		//check next chapter
		else if (typeof Trail.contents[+oldpos[0] + 1]['pages'] != 'undefined')
			newpos = (+oldpos[0] + 1) + '-' + 0;
		else {
			return;
		}
		return newpos;
	};

	var calcPrevPage = function trailPrevPage(Trail, oldpos) {
		if (oldpos == '0-0') return;

		oldpos = angular.isString(oldpos) ? oldpos.split("-") : [0, 1];
		var newpos;

		if (typeof Trail.contents[oldpos[0]]['pages'][+oldpos[1] - 1] != 'undefined')
			newpos = oldpos[0] + '-' + (+oldpos[1] - 1);
		else if (typeof Trail.contents[+oldpos[0] - 1]['pages'])
			newpos = (+oldpos[0] - 1) + '-' + (Trail.contents[+oldpos[0] - 1]['pages'].length - 1);
		else {
			return;
		}
		return newpos;
	};

	//go to the location of a trail page
	var activate = function trailActivate (indices) {
		//if passed nothing
		if (angular.isEmpty(indices)) {
			$location.search('position', null);
			return;
		}

		//if just pass chapter
		if (('' + indices).indexOf('-') < 0) {
			indices = indices + '-0';
		}

		//todo - check if already at url of trail

		$location.search('position', indices);
	};

	//icons for both page types and material types
	var trailIconMap = {
		'book'        : 'glyphicon glyphicon-book',
		'exercise'    : 'glyphicon glyphicon-edit',
		'eye'         : 'glyphicon glyphicon-eye-open',
		'info'        : 'glyphicon glyphicon-info-sign',
		'list'        : 'glyphicon glyphicon-list-alt',
		'picture'     : 'glyphicon glyphicon-picture',
		'quiz'        : 'glyphicon glyphicon-pencil',
		'schedule'    : 'glyphicon glyphicon-calendar',
		'slides'      : 'glyphicon glyphicon-th-large',
		'syllabus'    : 'glyphicon glyphicon-list-alt',
		'video'       : 'glyphicon glyphicon-film',
		'undefined'   : 'glyphicon glyphicon-file'         //fallthrough
	};

	var mapIcon = function trailMapIcon (iconName) {
		iconName = iconName || 'undefined';
		return trailIconMap[iconName]  || trailIconMap['undefined'];
	};

	var favorite = function trailFavorite (id) {
		console.log("favorite trail with id: " + id);
	};

	return {
		compile : compile,
		extractPage : extractPage,
		calcNextPage : calcNextPage,
		calcPrevPage : calcPrevPage,
		activate : activate,
		mapIcon : mapIcon,
		share : Clotho.share,
		favorite : favorite,
		downloadDependencies : downloadDependencies
	}
});
