'use strict';

angular.module('bioApp', [])
    .config(['$routeProvider', function ($routeProvider) {
        $routeProvider
            .when('/', {
                templateUrl:'bio/list-partial.html'
            })
            .when('/:id', {
                templateUrl:'bio/individual-partial.html',
                resolve : {
                    individual : function($http, $q, $route) {
                        var deferred = $q.defer();

                        $http({method: "GET", url : "bio/faculty/" + $route.current.params.id + ".json", cache: true })
                            .success(function (result) {
                                //testing - replace images while don't exist
                                result.thumb = placeholders.thumb;
                                result.headshot = placeholders.headshot;

                                deferred.resolve(result);
                            })
                            //note- couldn't find that person, view won't update
                            .error(function (result) {
                                deferred.reject("Person couldn't be found");
                            });

                        return deferred.promise;
                    }
                }
            })
            .otherwise({
                templateUrl:'bio/list-partial.html'
            })
    }])
    .run(['$rootScope', '$location', function ($rootScope, $location) {
        //note - reset the route if the individual isn't found...
        $rootScope.$on('$routeChangeError', function() {$location.path("/")})

    }])
    .controller('MenuCtrl', ['$scope', function($scope) {
        $scope.menuItems = [
            {"text" : "People", "url" : "/people"},
            {"text" : "Research", "url" : "/research"},
            {"text" : "Software", "url" : "/software"},
            {"text" : "Education", "url" : "/education"},
            {"text" : "Blog", "url" : "/blog"}
        ];

    }])
    .service('Faculty', ['$http', '$location', function($http, $location) {

        var getList = function() {
            //use array to maintain order, angular filters available e.g. in ng-repeat
            var list = {
                "head" : [
                        {"name" : "Chris Anderson", "id" : "chrisAnderson", "description" : "This is a short blurb about Chris Anderson and what he does in the lab. He is the boss of everything."},
                        {"name" : "Robert Frawley", "id" : "robertFrawley", "description" : "This is a short blurb about Robert and what he does...."}
                ],
                "present" : [
                    {
                        "title" : "Post Docs",
                        "members" :[
                            {"name":"David Sukovich","id":"davidSukovich"},
                            {"name":"Saurabh Srivastava","id":"saurabhSrivastava"}
                        ]
                    },
                    {
                        "title" : "Grad Students",
                        "members" :[
                            {"name":"Tim Hsiau","id":"timHsiau"},
                            {"name":"Gabriel Lopez","id":"gabrielLopez"}
                        ]
                    },
                    {
                        "title" : "Programmers",
                        "members" :[
                            {"name":"Max Bates","id":"maxBates"}
                        ]
                    },
                    {
                        "title" : "Undergrads",
                        "members" :[
                            {"name":"Paul Ruan","id":"paulRuan"},
                            {"name":"Jeff Tsui","id":"jeffTsui"},
                            {"name":"Jene Liang","id":"jeneLiang"}
                        ]
                    }
                ],
                "collab" : {
                    "title" : "Collaborators",
                    "filter" : "name",
                    "members" : [
                        {"name":"Douglas Densmore","id":"douglasDensmore"},
                        {"name":"Ras Bodik","id":"rasBodik"},
                        {"name":"Sanjit Seshia","id":"sanjitSeshia"},
                        {"name":"Sarah Chasins","id":"sarahChasins"},
                        {"name":"Ernst Oberortner","id":"ernstOberortner"},
                        {"name":"Stephanie Paige","id":"stephaniePaige"},
                        {"name":"Terry Johnson","id":"terryJohnson"},
                        {"name":"George Khushf","id":"georgeKhushf"}
                    ]
                },
                "past" : [
                    {
                        "title" : "Post Docs",
                        "members" :[
                            {"name":"Mariana Leguia","id":"marianaLeguia"},
                            {"name":"Phillip Elms","id":"phillipElms"}
                        ]
                    },
                    {
                        "title" : "Grad Students",
                        "members" :[
                            {"name":"Jin Huh","id":"jinHuh"},
                            {"name":"Josh Kittleson","id":"joshKittleson"}
                        ]
                    },
                    {
                        "title" : "Undergrads",
                        "members" :[
                            {"name":"Austin Day","id":"austinDay"},
                            {"name":"Neyer Christoph","id":"neyerChristoph"},
                            {"name":"David Tulga","id":"davidTulga"},
                            {"name":"Janna Serbo","id":"jannaSerbo"},
                            {"name":"Jenn Brophy","id":"jennBrophy"},
                            {"name":"Madhvi Venkatesh","id":"madhviVenkatesh"},
                            {"name":"Nima Emami","id":"nimaEmami"},
                            {"name":"Samantha Liang","id":"samanthaLiang"},
                            {"name":"Cheung Sherine","id":"cheungSherine"},
                            {"name":"Sushant Sundaresh","id":"sushantSundaresh"},
                            {"name":"Tahoura Samad","id":"tahouraSamad"},
                            {"name":"Hannah Cole","id":"hannahCole"},
                            {"name":"Shelly Cheng","id":"shellyCheng"}
                        ]
                    }
                ]
            };

            return list;

        };

        //generally, use logic in routing instead of this
        var getIndiv = function(id) {
            $http({method: "GET", url : "bio/faculty/" + id + ".json", cache: true });
        };

        var goHome = function() {
            $location.path("/");
        };

        var goToBio = function(id) {
            $location.path("/" + id);
        };

        return {
            getList: getList,
            getIndiv : getIndiv,
            goToBio : goToBio,
            goHome : goHome
        }
    }])
    .controller('ListCtrl', ['$scope', 'Faculty', function($scope, Faculty) {

        $scope.thumb = placeholders.thumb;
        $scope.large_thumb = placeholders.large_thumb;

        $scope.faculty = Faculty.getList();

        $scope.goToBio = Faculty.goToBio;

    }])
    .controller('IndividualCtrl', ['$scope', '$route', 'Faculty', function($scope, $route, Faculty) {

        //inherit from routeProvider.resolve(), promise fulfilled before render
        $scope.indiv = $route.current.locals.individual;

        $scope.id = $route.current.params.id;

        $scope.next = function() {

        };

        $scope.previous = function() {

        };

        $scope.home = function() {
            Faculty.goHome();
        }

    }]);


var placeholders = {};
placeholders.thumb = "data:image/png;base64,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";
placeholders.large_thumb = "data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAASwAAADICAYAAABS39xVAAAL60lEQVR4Xu3bCY8UVRcG4EI0ils0BIi7iQgiIhqNmrj9dnfiAlHQgAoqGOMScd8XvpxKqtP2NyNTyBvmwHMTA87UnD793MubW7drtpw9e/bcYBAgQKCBwBaB1WCWtEiAwCggsCwEAgTaCAisNlOlUQIEBJY1QIBAGwGB1WaqNEqAgMCyBggQaCMgsNpMlUYJEBBY1gABAm0EBFabqdIoAQICyxogQKCNgMBqM1UaJUBAYFkDBAi0ERBYbaZKowQICCxrgACBNgICq81UaZQAAYFlDRAg0EZAYLWZKo0SICCwrAECBNoICKw2U6VRAgQEljVAgEAbAYHVZqo0SoCAwLIGCBBoIyCw2kyVRgkQEFjWAAECbQQEVpup0igBAgLLGiBAoI2AwGozVRolQEBgWQMECLQREFhtpkqjBAgILGuAAIE2AgKrzVRplAABgWUNECDQRkBgtZkqjRIgILCsAQIE2ggIrDZTpVECBASWNUCAQBsBgdVmqjRKgIDAsgYIEGgjILDaTJVGCRAQWNYAAQJtBARWm6nSKAECAssaIECgjYDAajNVGiVAQGBZAwQItBEQWG2mSqMECAgsa4AAgTYCAqvNVGmUAAGBZQ0QINBGQGC1mSqNEiAgsKwBAgTaCAisNlOlUQIEBJY1QIBAGwGB1WaqNEqAgMCyBggQaCMgsNpMlUYJEBBY1gABAm0EBFabqdIoAQICyxogQKCNgMBqM1UaJUBAYFkDBAi0ERBYbaZKowQICCxrgACBNgICq81UaZQAAYFlDRAg0EZAYLWZKo0SICCwrAECBNoICKw2U6VRAgQEljVAgEAbAYHVZqo0SoCAwLIGCBBoIyCw2kyVRgkQEFjWAAECbQQEVpup0igBAgLLGiBAoI2AwGozVRolQEBgWQMECLQREFhtpkqjBAgILGuAAIE2AgKrzVRplAABgWUNECDQRkBgtZkqjRIgILCsAQIE2ggIrDZTpVECBASWNUCAQBsBgdVmqjRKgIDAsgYIEGgjILDaTJVGCRAQWNYAAQJtBARWm6nSKAECAssaIECgjYDAajNVGiVAQGBdRmvg77//Hv7666/hqquuGrZu3fqv72y6ti665pprLtq1czn/+OOP8UcuZg9z3tvcfl1/aQUE1qX1vyiv/ssvvwxHjx4dvv/++0W9a6+9dtizZ8+wa9euf7xGBcSxY8eGr7/+evH1q6++erjvvvuGu+6664KvnfNGfvvtt7GHb775ZvFjFVj333//cPvtt19wD3Pe25x+Xbt5BATW5pmLC+rk119/HV599dWhdhVrjQqBe+65Z/zWuXPnhpdffnmowFhrPPDAA8Odd945+9o5jf/555/DSy+9NO4E1xoVWA8++ODsHua8tzn9unZzCQiszTUfs7upndUXX3yx+LnaUX311VeLAKvbw+eff368Rfzss8+G9957b3Htzp07h++++24RYLXTevbZZ2dfO6fp48ePD59++uniR2655Zbhxx9/HCrIpvH4448P9fXN0O+c9+bavIDAyhvHXqF2Fa+99trw888/j6/x5JNPDjfddNNQt4i166rv13j66aeHbdu2DYcPHx7Onj07fm3ayfz+++/jrmvaoU1hsdFrr7vuuuHEiRPjuVntmm677bbFbWh9vXZzW7ZsGa6//vrxtvPQoUNjQNWoW9a77777/3Z++/fvH+tstIcKtznXxiZE4biAwIoTZ1/g7bffHn766afx0Pqxxx4bg6N2Ky+++OIYQhUWFVj1/boVm3YyU4hVd0eOHFmcadXtWO3SLvTaer3a0dX5VPU2jd27dw/33nvvWHe6JX3mmWeGCrwa77zzzvDll18uwnTv3r0X3MO/vbfVM7Ls7Kh+sQUE1sUWvYT1Kri+/fbb4cyZM4tdTAVC/QOuA+nlndRyWHz44YfDxx9/PHZeYVVhsdFrDxw4ML5W7ZymUbea9QFAna/VuPnmm4cnnnhi/Hv1WDuxCrYbb7xx/LN2gi+88MIiTGuHtX379g33MLffSzhFXvo/Cgis/wi4mX58+Xyodlq1w3r44YeHCpDlXVf1vBxY77///nD69OlFYO3bt2+xQzvftRVYNZZrLJtMO7xpJ7XqVUH65ptvjkE2jaeeemrceU27xPP1cCH9bqZ508vGBQTWxq02/ZVrhUYdpFc41UgGwHqf0i1/8rgKWAf+b7311j8+4axHK2rHlA7YTT+ZGlxTQGBdRgujdlT1D/3zzz8fdzzTON8t1nLQ1W6sQma9W8LVa2sHN4060K/D72lMt6O1y1odH3zwwfDJJ5/848t1AF8H8TVWPwxYb0f4X/q9jKb+inkrAqvxVNft1Lvvvju+gzpUrwPzKRyWPzWrc6m6bVp+/mk5AJYfjajntu64444NXzs941U91G5p+WHQ+lrdMq4+vLp8ZlbXVM+PPPLIeG41jdXntRL9Np76K7Z1gdV46msXUiFUt2N1ZvXcc88NdQtYox53mM6Fpk/o3njjjfG5qxrTA6WrB94HDx4cduzYMcy5turVQX89xrA6VvuqDwXqzGoat9566/Doo4+O/a+OOT3MubbxlF/xrQusxktgNWzq13Hqdq4ezFz+1ZsKhNq9fPTRR8PJkycX7/ihhx4aHzKdHjydHkmo0JtzbX0a+Morryye+6pnrqZnw+rFKgArCGssB0v9f71WhdX0zFj9WU/bV8jO6WHOtY2n/IpvXWA1XwKnTp0a6r/1Rj06UA+UVhid79dipgPvqrXRaytgKoSm32OcQq96mj55rHoVjhVcy893rddznUvV2dhGe5jTb/PpvuLbF1iXwRJY6wB72tnUGdLy7VY9Bf/666+Pz2Utjzq3qnOu5bGRa1cP2qfbz6pfu67pQdU6gK/grMP89X6PcHrt5V420sP0c3OuvQym/Yp8CwLrMpn2+oTwhx9+GIOodjk33HDD4inytd7idMtY4VG7sLqNW2/MuTbFOaeHOdem+lU3IyCwMq6qEiAQEBBYAVQlCRDICAisjKuqBAgEBARWAFVJAgQyAgIr46oqAQIBAYEVQFWSAIGMgMDKuKpKgEBAQGAFUJUkQCAjILAyrqoSIBAQEFgBVCUJEMgICKyMq6oECAQEBFYAVUkCBDICAivjqioBAgEBgRVAVZIAgYyAwMq4qkqAQEBAYAVQlSRAICMgsDKuqhIgEBAQWAFUJQkQyAgIrIyrqgQIBAQEVgBVSQIEMgICK+OqKgECAQGBFUBVkgCBjIDAyriqSoBAQEBgBVCVJEAgIyCwMq6qEiAQEBBYAVQlCRDICAisjKuqBAgEBARWAFVJAgQyAgIr46oqAQIBAYEVQFWSAIGMgMDKuKpKgEBAQGAFUJUkQCAjILAyrqoSIBAQEFgBVCUJEMgICKyMq6oECAQEBFYAVUkCBDICAivjqioBAgEBgRVAVZIAgYyAwMq4qkqAQEBAYAVQlSRAICMgsDKuqhIgEBAQWAFUJQkQyAgIrIyrqgQIBAQEVgBVSQIEMgICK+OqKgECAQGBFUBVkgCBjIDAyriqSoBAQEBgBVCVJEAgIyCwMq6qEiAQEBBYAVQlCRDICAisjKuqBAgEBARWAFVJAgQyAgIr46oqAQIBAYEVQFWSAIGMgMDKuKpKgEBAQGAFUJUkQCAjILAyrqoSIBAQEFgBVCUJEMgICKyMq6oECAQEBFYAVUkCBDICAivjqioBAgEBgRVAVZIAgYyAwMq4qkqAQEBAYAVQlSRAICMgsDKuqhIgEBAQWAFUJQkQyAgIrIyrqgQIBAQEVgBVSQIEMgICK+OqKgECAQGBFUBVkgCBjIDAyriqSoBAQEBgBVCVJEAgIyCwMq6qEiAQEBBYAVQlCRDICAisjKuqBAgEBARWAFVJAgQyAgIr46oqAQIBAYEVQFWSAIGMgMDKuKpKgEBAQGAFUJUkQCAjILAyrqoSIBAQEFgBVCUJEMgICKyMq6oECAQEBFYAVUkCBDICAivjqioBAgEBgRVAVZIAgYyAwMq4qkqAQEBAYAVQlSRAICMgsDKuqhIgEBAQWAFUJQkQyAgIrIyrqgQIBAQEVgBVSQIEMgL/A18IXjoHbg2UAAAAAElFTkSuQmCC";
placeholders.headshot = "data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAWgAAAEOCAYAAACkSI2SAAARg0lEQVR4Xu3chY8kVRAH4D5cg7uEENzd5W/HHYK7a5DDXVOd9KSvmZXb45aqq68TAuxOT9f76uU3vW/ezJ69e/f+PTgIECBAIJ3AHgGdricKIkCAwCggoE0EAgQIJBUQ0EkboywCBAgIaHOAAAECSQUEdNLGKIsAAQIC2hwgQIBAUgEBnbQxyiJAgICANgcIECCQVEBAJ22MsggQICCgzQECBAgkFRDQSRujLAIECAhoc4AAAQJJBQR00sYoiwABAgLaHCBAgEBSAQGdtDHKIkCAgIA2BwgQIJBUQEAnbYyyCBAgIKDNAQIECCQVENBJG6MsAgQICGhzgAABAkkFBHTSxiiLAAECAtocIECAQFIBAZ20McoiQICAgDYHCBAgkFRAQCdtjLIIECAgoM0BAgQIJBUQ0EkboywCBAgIaHOAAAECSQUEdNLGKIsAAQIC2hwgQIBAUgEBnbQxyiJAgICANgcIECCQVEBAJ22MsggQICCgzQECBAgkFRDQSRujLAIECAhoc4AAAQJJBQR00sYoiwABAgLaHCBAgEBSAQGdtDHKIkCAgIA2BwgQIJBUQEAnbYyyCBAgIKDNAQIECCQVENBJG6MsAgQICGhzgAABAkkFBHTSxiiLAAECAtocIECAQFIBAZ20McoiQICAgDYHCBAgkFRAQCdtjLIIECAgoM0BAgQIJBUQ0EkboywCBAgIaHOAAAECSQUEdNLGKIsAAQIC2hwgQIBAUgEBnbQxyiJAgICANgcIECCQVEBAJ22MsggQICCgzQECBAgkFRDQSRujLAIECAhoc4AAAQJJBQR00sYoiwABAgLaHCBAgEBSAQGdtDHKIkCAgIA2BwgQIJBUQEAnbYyyCBAgIKDNAQIECCQVENBJG6MsAgQICGhzgAABAkkFBHTSxiiLAAECAtocIECAQFIBAZ20McoiQICAgDYHCBAgkFRAQCdtjLIIECAgoM0BAgQIJBUQ0EkboywCBAgIaHOAAAECSQUEdNLGKIsAAQIC2hwgQIBAUgEBnbQxyiJAgICANgcIECCQVEBAJ22MsggQICCgzQECBAgkFRDQSRujLAIECAhoc4AAAQJJBQR00sYoiwABAgLaHCBAgEBSAQGdtDHKIkCAgIA2BwgQIJBUQEAnbYyyCBAgIKDNAQIECCQVENBJG6MsAgQICGhzgAABAkkFBHTSxiiLAAECAtocIECAQFIBAZ20McoiQICAgDYHCBAgkFRAQCdtjLIIECAgoM0BAgQIJBUQ0EkboywCBAgIaHOAAAECSQUEdNLGKIsAAQIC2hwgQIBAUgEBnbQxyiJAgICANgcIECCQVEBAJ22MsggQICCgzQECBAgkFRDQSRujLAIECAhoc4AAAQJJBQR00sYoiwABAgLaHCBAgEBSAQGdtDHKIkCAgIA2BwgQIJBUQEAnbYyyCBAgIKDNAQIECCQVENBJG6MsAgQICGhzgAABAkkFBHTSxiiLAAECAtocIECAQFIBAZ20McoiQICAgDYHCBAgkFRAQCdtjLIIECAgoM0BAgQIJBUQ0EkboywCBAgIaHOAAAECSQUEdNLGKIsAAQIC2hzYUuC7774b3n777SH+/ddff42PP/roo4cLLrhg/Gej45133hk++OCD4Y8//hgfsmfPnuGMM84YLr/88vH85fHzzz8Pb7311vDll18Of/755/jrY445ZjjrrLOGSy65ZDz/vzh++OGH4d133x2+/vrr8ToxpqjnvPPOGy666KLVdaLuV155ZduXjOe47LLLVufv1ni2XaAHlhMQ0OVatrsFf/zxx8Orr7664UWPPPLI4e677x6OOOKI1WP+/vvv4dFHHx0ioNYdEbQ33njjcOqpp65+HeH/1FNPDXHuuiOuc++99w6HHXbYAQF88cUXw/PPP7/hcxx++OHD7bffPhx33HHD77//Pjz00EOrF6WtLhwG999//xjQuzWerWry+9oCArp2/w5q9RGwEbQbheZ08Qjam266aVVLBGAE4WZHhNl99903Bm48fwThb7/9tuk5Z5555nDdddfteMy//vrr8Mgjj2wZuEcdddT4YhB30PsT0NN5UeBujGfHEE4sIyCgy7Rq9wt98803h/fff3914dNOO21cnoi76vnP447xnnvuGZcJYvng8ccfX50z3S3HHfDTTz+9WrqIB9xwww3D6aefPuzdu3d49tln9znntttuG3788cfhpZde2ufnEerxXDs5YlkjlmqmI+6Wr7322vE6sbQyfyGKn8dyzDPPPLM20OOF5dtvv92njBNPPHGIumPpZDfGsxMD59QSENC1+rWr1UbQRuDGMf/zPf4/giuCKI4I4bvuums49thj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";